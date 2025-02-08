// The main ray tracer.

#pragma warning(disable : 4786)

#include "RayTracer.h"
#include "scene/light.h"
#include "scene/material.h"
#include "scene/ray.h"

#include "parser/JsonParser.h"
#include "parser/Parser.h"
#include "parser/Tokenizer.h"
#include <json.hpp>

#include "ui/TraceUI.h"
#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>
#include <stdbool.h>
#include <glm/gtx/io.hpp>
#include <string.h> // for memset

#include <fstream>
#include <iostream>

using namespace std;
extern TraceUI *traceUI;

// Use this variable to decide if you want to print out debugging messages. Gets
// set in the "trace single ray" mode in TraceGLWindow, for example.
bool debugMode = false;

// Trace a top-level ray through pixel(i,j), i.e. normalized window coordinates
// (x,y), through the projection plane, and out into the scene. All we do is
// enter the main ray-tracing method, getting things started by plugging in an
// initial ray weight of (0.0,0.0,0.0) and an initial recursion depth of 0.

glm::dvec3 RayTracer::trace(double x, double y) {
  // Clear out the ray cache in the scene for debugging purposes,
  if (TraceUI::m_debug) {
    scene->clearIntersectCache();
  }

  ray r(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0), glm::dvec3(1, 1, 1),
        ray::VISIBILITY);
  scene->getCamera().rayThrough(x, y, r);
  double dummy;
  glm::dvec3 ret =
      traceRay(r, glm::dvec3(1.0, 1.0, 1.0), traceUI->getDepth(), dummy);
  ret = glm::clamp(ret, 0.0, 1.0);
  return ret;
}

glm::dvec3 RayTracer::tracePixel(int i, int j) {
  glm::dvec3 col(0, 0, 0);

  if (!sceneLoaded())
    return col;

  double x = double(i) / double(buffer_width);
  double y = double(j) / double(buffer_height);

  unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;
  col = trace(x, y);

  pixel[0] = (int)(255.0 * col[0]);
  pixel[1] = (int)(255.0 * col[1]);
  pixel[2] = (int)(255.0 * col[2]);
  return col;
}

#define VERBOSE 0

// Do recursive ray tracing! You'll want to insert a lot of code here (or places
// called from here) to handle reflection, refraction, etc etc.
glm::dvec3 RayTracer::traceRay(ray &r, const glm::dvec3 &thresh, int depth,
                               double &t) {
  isect i;
  glm::dvec3 colorC;
#if VERBOSE
  std::cerr << "== current depth: " << depth << std::endl;
#endif

  if (scene->intersect(r, i)) {
    // YOUR CODE HERE

    // An intersection occurred!  We've got work to do. For now, this code gets
    // the material for the surface that was intersected, and asks that material
    // to provide a color for the ray.

    // This is a great place to insert code for recursive ray tracing. Instead
    // of just returning the result of shade(), add some more steps: add in the
    // contributions from reflected and refracted rays.

    // An intersection occurred.
    const Material &m = i.getMaterial();
    
    // Start with the local (base) shading contribution.
    colorC = m.shade(scene.get(), r, i);

    // If we have remaining recursion depth, add reflective and refractive contributions.
    if (depth > 0) {
        // --- Reflection ---
        // Get the normalized surface normal at the intersection.
        glm::dvec3 N = glm::normalize(i.getN());
        // Get the normalized incident direction (from the ray).
        glm::dvec3 I = glm::normalize(r.getDirection());
        // Compute the reflection direction: R = I - 2 (I·N) N.
        glm::dvec3 reflectDir = I - 2.0 * glm::dot(I, N) * N;
        reflectDir = glm::normalize(reflectDir);
        
        // Offset the intersection point slightly along the normal to avoid self-intersection.
        glm::dvec3 reflectOrigin = r.at(i) + RAY_EPSILON * N;
        // Construct a reflection ray.
        ray reflectRay(reflectOrigin, reflectDir, r.getAtten(), ray::REFLECTION);
        
        // Recursively trace the reflection ray.
        double tReflect;
        glm::dvec3 reflectColor = traceRay(reflectRay, thresh, depth - 1, tReflect);
        
        // Get the reflectivity coefficient from the material (assumed to be in [0,1]).
        // Add the reflection contribution.
        colorC += m.kr(i) * reflectColor;

        // --- Refraction ---
        // Check if the material is transparent.
        if (m.Trans() && depth > 0) {
        // Get the material's index of refraction.
        double materialIndex = m.index(i);
        
        // Get normalized incident direction and surface normal.
        glm::dvec3 I = glm::normalize(r.getDirection());
        glm::dvec3 N = glm::normalize(i.getN());
        
        // Determine if the ray is entering or exiting the material.
        // When entering, eta = air_index / material_index, and when exiting, eta = material_index / air_index.
        // (We assume the surrounding medium is air with index 1.0.)
        double eta;
        glm::dvec3 effectiveN = N; // effective normal used for refraction calculation.
        if (glm::dot(I, N) < 0) {
            // Ray is entering the medium.
            eta = 1.0 / materialIndex;
        } else {
            // Ray is leaving the medium; flip the normal.
            effectiveN = -N;
            eta = materialIndex;
        }
        
        // Compute the refracted direction using GLM's built-in function.
        glm::dvec3 refractDir = glm::refract(I, effectiveN, eta);
        
        // Check for total internal reflection. glm::refract returns a zero vector when no refraction occurs.
        if (glm::length(refractDir) > 1e-6) {
            refractDir = glm::normalize(refractDir);
            // Offset the ray origin slightly to avoid self-intersections.
            glm::dvec3 refractOrigin = r.at(i) - RAY_EPSILON * N;
            // Construct the refraction ray.
            ray refractRay(refractOrigin, refractDir, r.getAtten(), ray::REFRACTION);
            double tRefract;
            glm::dvec3 refractColor = traceRay(refractRay, thresh, depth - 1, tRefract);
            
            // Add the refraction contribution to the final color.
            colorC += m.Trans() * refractColor;
        }
      }
    }
  }   // const Material &m = i.getMaterial();
    // colorC = m.shade(scene.get(), r, i);
  else {
    // No intersection. This ray travels to infinity, so we color
    // it according to the background color, which in this (simple)
    // case is just black.
    //
    // FIXME: Add CubeMap support here.
    // TIPS: CubeMap object can be fetched from
    // traceUI->getCubeMap();
    //       Check traceUI->cubeMap() to see if cubeMap is loaded
    //       and enabled.

    // CubeMap *cube = traceUI->cubeMap();
    // if (cube /* && cube->isEnabled() */) { // Uncomment & implement isEnabled() if desired.
    //     colorC = cube->getColor(r);
    // } else {
        // Fallback: return a default background color (e.g., light blue).
        
        if(traceUI->cubeMap()){
        CubeMap* cube_map = traceUI->getCubeMap();
        colorC = cube_map->getColor(r);			
      } else{
        colorC = glm::dvec3(0, 0, 0);
      }
    }
  
#if VERBOSE
  std::cerr << "== depth: " << depth + 1 << " done, returning: " << colorC
            << std::endl;
#endif
  return colorC;
}

RayTracer::RayTracer()
    : scene(nullptr), buffer(0), thresh(0), buffer_width(0), buffer_height(0),
      m_bBufferReady(false) {
}

RayTracer::~RayTracer() {}

void RayTracer::getBuffer(unsigned char *&buf, int &w, int &h) {
  buf = buffer.data();
  w = buffer_width;
  h = buffer_height;
}

double RayTracer::aspectRatio() {
  return sceneLoaded() ? scene->getCamera().getAspectRatio() : 1;
}

bool RayTracer::loadScene(const char *fn) {
  ifstream ifs(fn);
  if (!ifs) {
    string msg("Error: couldn't read scene file ");
    msg.append(fn);
    traceUI->alert(msg);
    return false;
  }

  // Check if fn ends in '.ray'
  bool isRay = false;
  const char *ext = strrchr(fn, '.');
  if (ext && !strcmp(ext, ".ray"))
    isRay = true;

  // Strip off filename, leaving only the path:
  string path(fn);
  if (path.find_last_of("\\/") == string::npos)
    path = ".";
  else
    path = path.substr(0, path.find_last_of("\\/"));

  if (isRay) {
    // .ray Parsing Path
    // Call this with 'true' for debug output from the tokenizer
    Tokenizer tokenizer(ifs, false);
    Parser parser(tokenizer, path);
    try {
      scene.reset(parser.parseScene());
    } catch (SyntaxErrorException &pe) {
      traceUI->alert(pe.formattedMessage());
      return false;
    } catch (ParserException &pe) {
      string msg("Parser: fatal exception ");
      msg.append(pe.message());
      traceUI->alert(msg);
      return false;
    } catch (TextureMapException e) {
      string msg("Texture mapping exception: ");
      msg.append(e.message());
      traceUI->alert(msg);
      return false;
    }
  } else {
    // JSON Parsing Path
    try {
      JsonParser parser(path, ifs);
      scene.reset(parser.parseScene());
    } catch (ParserException &pe) {
      string msg("Parser: fatal exception ");
      msg.append(pe.message());
      traceUI->alert(msg);
      return false;
    } catch (const json::exception &je) {
      string msg("Invalid JSON encountered ");
      msg.append(je.what());
      traceUI->alert(msg);
      return false;
    }
  }

  if (!sceneLoaded())
    return false;

  return true;
}

bool RayTracer::refract(const glm::dvec3 &I, const glm::dvec3 &N, double n, glm::dvec3 &T) const {
    // I: incident direction (normalized)
    // N: surface normal at the intersection (normalized, pointing outwards)
    // n: index of refraction of the material (e.g., 1.5)
    // T: computed refracted direction (if refraction occurs)
    //
    // We assume that the external medium (air) has index 1.0.
    
    // Compute the cosine of the angle between I and N.
    double cosi = glm::dot(I, N);
    double etai = 1.0, etat = n;
    glm::dvec3 n_normal = N;  // This will be the effective normal.

    // If cosi is positive, the ray is exiting the material.
    if (cosi > 0) {
        // Swap indices and flip the normal.
        std::swap(etai, etat);
        n_normal = -N;
    } else {
        // If the ray is entering, make cosi positive.
        cosi = -cosi;
    }
    
    // Compute the ratio of indices.
    double eta = etai / etat;
    
    // Compute k = 1 - eta^2*(1 - cosi^2)
    double k = 1.0 - eta * eta * (1.0 - cosi * cosi);
    if (k < 0.0) {
        // Total internal reflection; no refraction occurs.
        return false;
    }
    
    // Compute the refracted direction.
    T = eta * I + (eta * cosi - sqrt(k)) * n_normal;
    T = glm::normalize(T);  // Ensure normalization.
    return true;
}

void RayTracer::traceSetup(int w, int h) {
  size_t newBufferSize = w * h * 3;
  if (newBufferSize != buffer.size()) {
    bufferSize = newBufferSize;
    buffer.resize(bufferSize);
  }
  buffer_width = w;
  buffer_height = h;
  std::fill(buffer.begin(), buffer.end(), 0);
  m_bBufferReady = true;

  /*
   * Sync with TraceUI
   */

  threads = traceUI->getThreads();
  block_size = traceUI->getBlockSize();
  thresh = traceUI->getThreshold();
  samples = traceUI->getSuperSamples();
  aaThresh = traceUI->getAaThreshold();

  // YOUR CODE HERE
  // FIXME: Additional initializations
}

/*
 * RayTracer::traceImage
 *
 *	Trace the image and store the pixel data in RayTracer::buffer.
 *
 *	Arguments:
 *		w:	width of the image buffer
 *		h:	height of the image buffer
 *
 */
void RayTracer::traceImage(int w, int h) {
  // Always call traceSetup before rendering anything.
  traceSetup(w, h);
  for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			tracePixel(i, j);
		}
	}
  // YOUR CODE HERE
  // FIXME: Start one or more threads for ray tracing
  //
  // TIPS: Ideally, the traceImage should be executed asynchronously,
  //       i.e. returns IMMEDIATELY after working threads are launched.
  //
  //       An asynchronous traceImage lets the GUI update your results
  //       while rendering.
}

int RayTracer::aaImage() {
  // YOUR CODE HERE
  // FIXME: Implement Anti-aliasing here
  //
  // TIP: samples and aaThresh have been synchronized with TraceUI by
  //      RayTracer::traceSetup() function
  // Loop over every pixel in the image.
    for (int j = 0; j < buffer_height; ++j) {
        for (int i = 0; i < buffer_width; ++i) {
            // Compute normalized coordinate bounds for the pixel.
            // (Assumes x and y in [0,1] across the image.)
            double x0 = static_cast<double>(i) / buffer_width;
            double y0 = static_cast<double>(j) / buffer_height;
            double x1 = static_cast<double>(i + 1) / buffer_width;
            double y1 = static_cast<double>(j + 1) / buffer_height;
            
            // Use the 'samples' parameter (synchronized with TraceUI) as the maximum recursion depth.
            int maxDepth = samples;
            
            // Compute the anti–aliased color for this pixel.
            glm::dvec3 color = adaptiveSample(x0, y0, x1, y1, maxDepth);
            color = glm::clamp(color, 0.0, 1.0);
            
            // Write the computed color into the buffer.
            setPixel(i, j, color);
        }
    }
    
    // Return a status code (e.g., 1 indicates success).
    return 1;
}

// Helper function: recursively sample a pixel region defined by the normalized coordinates 
// [x0, x1] x [y0, y1].
// 'depth' is the remaining recursion depth (typically initialized from 'samples').
// 'aaThresh' is the minimum color difference required to trigger further subdivision.
glm::dvec3 RayTracer::adaptiveSample(double x0, double y0, double x1, double y1, int depth) {
    // Sample the four corners.
    glm::dvec3 c00 = trace(x0, y0);
    glm::dvec3 c10 = trace(x1, y0);
    glm::dvec3 c01 = trace(x0, y1);
    glm::dvec3 c11 = trace(x1, y1);
    // Sample the center.
    glm::dvec3 center = trace((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    
    // Compute the average of the corner colors.
    glm::dvec3 avgCorners = (c00 + c10 + c01 + c11) / 4.0;
    
    // If we've reached the maximum recursion depth or the color difference is small,
    // return the average of the five samples (four corners and the center).
    if (depth <= 0 || glm::length(center - avgCorners) < aaThresh) {
        return (c00 + c10 + c01 + c11 + center) / 5.0;
    }
    
    // Otherwise, subdivide the region into four quadrants.
    double xm = (x0 + x1) / 2.0;
    double ym = (y0 + y1) / 2.0;
    glm::dvec3 s1 = adaptiveSample(x0, y0, xm, ym, depth - 1);
    glm::dvec3 s2 = adaptiveSample(xm, y0, x1, ym, depth - 1);
    glm::dvec3 s3 = adaptiveSample(x0, ym, xm, y1, depth - 1);
    glm::dvec3 s4 = adaptiveSample(xm, ym, x1, y1, depth - 1);
    
    return (s1 + s2 + s3 + s4) / 4.0;
}

bool RayTracer::checkRender() {
  // YOUR CODE HERE
  // FIXME: Return true if tracing is done.
  //        This is a helper routine for GUI.
  //
  // TIPS: Introduce an array to track the status of each worker thread.
  //       This array is maintained by the worker threads.
  return true;
}

void RayTracer::waitRender() {
  // YOUR CODE HERE
  // FIXME: Wait until the rendering process is done.
  //        This function is essential if you are using an asynchronous
  //        traceImage implementation.
  //
  // TIPS: Join all worker threads here.
}


glm::dvec3 RayTracer::getPixel(int i, int j) {
  unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;
  return glm::dvec3((double)pixel[0] / 255.0, (double)pixel[1] / 255.0,
                    (double)pixel[2] / 255.0);
}

void RayTracer::setPixel(int i, int j, glm::dvec3 color) {
  unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;

  pixel[0] = (int)(255.0 * color[0]);
  pixel[1] = (int)(255.0 * color[1]);
  pixel[2] = (int)(255.0 * color[2]);
}
