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

glm::dvec3 RayTracer::trace(double x, double y)
{
  // Clear out the ray cache in the scene for debugging purposes,
  if (TraceUI::m_debug)
  {
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

glm::dvec3 RayTracer::tracePixel(int i, int j)
{
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
                               double &t)
{
  isect i;
  glm::dvec3 colorC;
#if VERBOSE
  std::cerr << "== current depth: " << depth << std::endl;
#endif

  if (scene->intersect(r, i))
  {
    // YOUR CODE HERE

    // An intersection occurred!  We've got work to do. For now, this code gets
    // the material for the surface that was intersected, and asks that material
    // to provide a color for the ray.

    // This is a great place to insert code for recursive ray tracing. Instead
    // of just returning the result of shade(), add some more steps: add in the
    // contributions from reflected and refracted rays.

    // An intersection occurred.
    //--------jerry code-------
    const Material &m = i.getMaterial();

    // Start with the local (base) shading contribution.
    colorC = m.shade(scene.get(), r, i);

    // If we have remaining recursion depth, add reflective and refractive contributions.
    if (depth > 0)
    {
      // --- Reflection ---
      // Get the normalized surface normal at the intersection.
      glm::dvec3 N = glm::normalize(i.getN());
      // Get the normalized incident direction (from the ray).
      glm::dvec3 I = glm::normalize(r.getDirection());
      // Compute the reflection direction: R = I - 2 (I·N) N.
      glm::dvec3 reflectDir = I - 2.0 * glm::dot(I, N) * N;
      reflectDir = glm::normalize(reflectDir);

      // Offset the intersection point slightly along the normal to avoid self-intersection.
      glm::dvec3 reflectOrigin = r.at(i);
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

      if (m.Trans() && depth > 0)
      {
        // Get the material's index of refraction.
        double materialIndex = m.index(i);

        // Get normalized incident direction and surface normal.
        glm::dvec3 I = glm::normalize(r.getDirection());
        glm::dvec3 N = glm::normalize(i.getN());

        // Determine if the ray is entering or exiting the material.
        // When entering, use eta = air_index / materialIndex;
        // When exiting, flip the normal and use eta = materialIndex / air_index.
        double eta;
        glm::dvec3 effectiveN = N; // effective normal for refraction computation
        bool entering = true;
        if (glm::dot(I, N) <= 0.0)
        {
          // Ray is entering the medium.
          eta = 1.0 / materialIndex;
          entering = true;
        }
        else
        {
          // Ray is exiting the medium; flip the normal.
          effectiveN = -N;
          eta = materialIndex; // equivalent to materialIndex / 1.0
          entering = false;
        }
        // double refracSquared = 1 - materialIndex * materialIndex * (1 - glm::dot(I, N) * glm::dot(I, N));

        // If the ray is nearly parallel or the material is nearly air, just continue in the same direction.
        // if (fabs(glm::dot(I, N)) >= 0.999 || fabs(materialIndex - 1.0) < 1e-6)
        // {
        //   ray rRefract(r.at(i) + (entering ? -RAY_EPSILON : RAY_EPSILON) * N, I, r.getAtten(), ray::REFRACTION);
        //   double tRefract;
        //   glm::dvec3 refractColor = traceRay(rRefract, thresh, depth - 1, tRefract);
        //   colorC += m.kt(i) * refractColor;
        // }
        // if (refracSquared > 0)
        // {
        // Compute the refracted direction.
        glm::dvec3 refractDir = glm::refract(I, effectiveN, eta);

        // Check for total internal reflection.
        if (glm::length(refractDir) > 1e-6)
        {
          refractDir = glm::normalize(refractDir);

          // Offset the ray origin to avoid self-intersection.
          // When entering, offset *into* the object;
          // when exiting, offset outward.
          glm::dvec3 refractOrigin;
          if (entering)
            refractOrigin = r.at(i) - RAY_EPSILON * N;
          else
            refractOrigin = r.at(i) + RAY_EPSILON * N;

          ray refractRay(refractOrigin, refractDir, r.getAtten(), ray::REFRACTION);

          // Optionally, you can compute an attenuation based on how far the ray travels inside the medium.
          double distance = 0.0;
          isect i2;
          if (scene->intersect(refractRay, i2))
          {
            distance = i2.getT();
          }
          // glm::dvec3 attenuation = glm::pow(m.kt(i), glm::dvec3(distance));
          glm::dvec3 attenuation = entering ? glm::pow(m.kt(i), glm::dvec3(distance)) : glm::dvec3(1.0, 1.0, 1.0);

          // if (distance == 0.0)
          // {
          //   printf("attenuation: %f\n", attenuation);
          // }
          double tRefract;
          glm::dvec3 refractColor = traceRay(refractRay, thresh, depth - 1, tRefract);

          // Combine the refracted contribution.
          colorC += attenuation * refractColor;
        }
        // }
        // else
        // {
        //   return colorC;
        // }
      }
    }
  } // const Material &m = i.getMaterial();
    // colorC = m.shade(scene.get(), r, i);
  else
  {
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

    if (traceUI->cubeMap())
    {
      CubeMap *cube_map = traceUI->getCubeMap();
      colorC = cube_map->getColor(r);
    }
    else
    {
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
      m_bBufferReady(false)
{
}

RayTracer::~RayTracer() {}

void RayTracer::getBuffer(unsigned char *&buf, int &w, int &h)
{
  buf = buffer.data();
  w = buffer_width;
  h = buffer_height;
}

double RayTracer::aspectRatio()
{
  return sceneLoaded() ? scene->getCamera().getAspectRatio() : 1;
}

bool RayTracer::loadScene(const char *fn)
{
  ifstream ifs(fn);
  if (!ifs)
  {
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

  if (isRay)
  {
    // .ray Parsing Path
    // Call this with 'true' for debug output from the tokenizer
    Tokenizer tokenizer(ifs, false);
    Parser parser(tokenizer, path);
    try
    {
      scene.reset(parser.parseScene());
    }
    catch (SyntaxErrorException &pe)
    {
      traceUI->alert(pe.formattedMessage());
      return false;
    }
    catch (ParserException &pe)
    {
      string msg("Parser: fatal exception ");
      msg.append(pe.message());
      traceUI->alert(msg);
      return false;
    }
    catch (TextureMapException e)
    {
      string msg("Texture mapping exception: ");
      msg.append(e.message());
      traceUI->alert(msg);
      return false;
    }
  }
  else
  {
    // JSON Parsing Path
    try
    {
      JsonParser parser(path, ifs);
      scene.reset(parser.parseScene());
    }
    catch (ParserException &pe)
    {
      string msg("Parser: fatal exception ");
      msg.append(pe.message());
      traceUI->alert(msg);
      return false;
    }
    catch (const json::exception &je)
    {
      string msg("Invalid JSON encountered ");
      msg.append(je.what());
      traceUI->alert(msg);
      return false;
    }
  }

  if (!sceneLoaded())
    return false;

  scene->buildKdTree();
  
  return true;
}

bool RayTracer::refract(const glm::dvec3 &I, const glm::dvec3 &N, double n, glm::dvec3 &T) const
{
  // I: incident direction (normalized)
  // N: surface normal at the intersection (normalized, pointing outwards)
  // n: index of refraction of the material (e.g., 1.5)
  // T: computed refracted direction (if refraction occurs)
  //
  // We assume that the external medium (air) has index 1.0.

  // Compute the cosine of the angle between I and N.
  double cosi = glm::dot(I, N);
  double etai = 1.0, etat = n;
  glm::dvec3 n_normal = N; // This will be the effective normal.

  // If cosi is positive, the ray is exiting the material.
  if (cosi > 0)
  {
    // Swap indices and flip the normal.
    std::swap(etai, etat);
    n_normal = -N;
  }
  else
  {
    // If the ray is entering, make cosi positive.
    cosi = -cosi;
  }

  // Compute the ratio of indices.
  double eta = etai / etat;

  // Compute k = 1 - eta^2*(1 - cosi^2)
  double k = 1.0 - eta * eta * (1.0 - cosi * cosi);
  if (k < 0.0)
  {
    // Total internal reflection; no refraction occurs.
    return false;
  }

  // Compute the refracted direction.
  T = eta * I + (eta * cosi - sqrt(k)) * n_normal;
  T = glm::normalize(T); // Ensure normalization.
  return true;
}

void RayTracer::traceSetup(int w, int h)
{
  size_t newBufferSize = w * h * 3;
  if (newBufferSize != buffer.size())
  {
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
void RayTracer::traceImage(int w, int h)
{
  // Always call traceSetup before rendering anything.
  traceSetup(w, h);
  for (int i = 0; i < w; i++)
  {
    for (int j = 0; j < h; j++)
    {
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

int RayTracer::aaImage()
{
  // YOUR CODE HERE
  // FIXME: Implement Anti-aliasing here
  //
  // TIP: samples and aaThresh have been synchronized with TraceUI by
  //      RayTracer::traceSetup() function
  // Loop over every pixel in the image.
  // If only one sample is specified, no anti-aliasing is needed.
  if (samples <= 1)
    return 1;

  // Loop over each output pixel.
  for (int j = 0; j < buffer_height; j++)
  {
    for (int i = 0; i < buffer_width; i++)
    {
      glm::dvec3 accumColor(0, 0, 0);

      // For each pixel, subdivide it into a grid of sub-samples.
      // 'samples' is the number of sub-samples per dimension.
      for (int sy = 0; sy < samples; sy++)
      {
        for (int sx = 0; sx < samples; sx++)
        {
          // Compute normalized coordinates for the sub-sample.
          // Each pixel covers [i/width, (i+1)/width] horizontally,
          // and [j/height, (j+1)/height] vertically.
          // We choose the center of each subpixel.
          double subX = (i + (sx + 0.5) / samples) / buffer_width;
          double subY = (j + (sy + 0.5) / samples) / buffer_height;

          // Trace a ray through the subpixel and accumulate its color.
          accumColor += trace(subX, subY);
        }
      }

      // Average the colors from all sub-samples.
      accumColor /= (samples * samples);

      // Write the computed color into the buffer.
      unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;
      pixel[0] = (int)(255.0 * glm::clamp(accumColor[0], 0.0, 1.0));
      pixel[1] = (int)(255.0 * glm::clamp(accumColor[1], 0.0, 1.0));
      pixel[2] = (int)(255.0 * glm::clamp(accumColor[2], 0.0, 1.0));
    }
  }

  return 1;
}

bool RayTracer::checkRender()
{
  // YOUR CODE HERE
  // FIXME: Return true if tracing is done.
  //        This is a helper routine for GUI.
  //
  // TIPS: Introduce an array to track the status of each worker thread.
  //       This array is maintained by the worker threads.
  return true;
}

void RayTracer::waitRender()
{
  // YOUR CODE HERE
  // FIXME: Wait until the rendering process is done.
  //        This function is essential if you are using an asynchronous
  //        traceImage implementation.
  //
  // TIPS: Join all worker threads here.
}

glm::dvec3 RayTracer::getPixel(int i, int j)
{
  unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;
  return glm::dvec3((double)pixel[0] / 255.0, (double)pixel[1] / 255.0,
                    (double)pixel[2] / 255.0);
}

void RayTracer::setPixel(int i, int j, glm::dvec3 color)
{
  unsigned char *pixel = buffer.data() + (i + j * buffer_width) * 3;

  pixel[0] = (int)(255.0 * color[0]);
  pixel[1] = (int)(255.0 * color[1]);
  pixel[2] = (int)(255.0 * color[2]);
}