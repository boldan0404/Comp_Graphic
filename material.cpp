#include "material.h"
#include "../ui/TraceUI.h"
#include "light.h"
#include "ray.h"
extern TraceUI *traceUI;

#include "../fileio/images.h"
#include <glm/gtx/io.hpp>
#include <iostream>

using namespace std;
extern bool debugMode;

Material::~Material() {}

// Apply the phong model to this point on the surface of the object, returning
// the color of that point.
glm::dvec3 Material::shade(Scene *scene, const ray &r, const isect &i) const
{
  // YOUR CODE HERE

  // For now, this method just returns the diffuse color of the object.
  // This gives a single matte color for every distinct surface in the
  // scene, and that's it.  Simple, but enough to get you started.
  // (It's also inconsistent with the phong model...)

  // Your mission is to fill in this method with the rest of the phong
  // shading model, including the contributions of all the light sources.
  // You will need to call both distanceAttenuation() and
  // shadowAttenuation()
  // somewhere in your code in order to compute shadows and light falloff.
  //	if( debugMode )
  //		std::cout << "Debugging Phong code..." << std::endl;

  // When you're iterating through the lights,
  // you'll want to use code that looks something
  // like this:
  //
  // for ( const auto& pLight : scene->getAllLights() )
  // {
  //              // pLight has type Light*
  // 		.
  // 		.
  // 		.
  // }
  // return kd(i);

  // different color channels for that equation, different light source same equation?
  // giving all the light sources, what will be the color at the interction point

  // initialize the color
  // glm::dvec3 colorC(0.0, 0.0, 0.0);
  // glm::dvec3 kd = kd(i);
  // glm::dvec3 ks = ks(i);
  // double alpha = shininess(i);

  // glm::dvec3 ambientTerm = ka(i) * scene->ambient();

  // // intersection point
  // glm::dvec3 p = r.at(i);

  // // add the ambient term and emissive term
  // colorC += ambientTerm + ke(i);

  // for (const auto &pLight : scene->getAllLights())
  // {
  //   // diffuse term:
  //   glm::dvec3 n = glm::normalize(i.getN());
  //   glm::dvec3 l = glm::normalize(pLight->getDirection(p)); // l: from intersection point to light source
  //   double lDotN = glm::max(glm::dot(l, n), 0.0);
  //   glm::dvec3 diffuseTerm = kd * lDotN * pLight->getColor();

  //   // specular term:
  //   glm::dvec3 r = glm::normalize(-l + 2 * glm::dot(-l, n) * n);
  //   glm::dvec3 v = glm::normalize(-r.getDirection());
  //   double vDotR = glm::max(glm::dot(r, v), 0.0);
  //   glm::dvec3 specularTerm = ks * pow(vDotR, alpha) * pLight->getColor();

  //   // Attenuation
  //   double distanceAtten = pLight->distanceAttenuation(p);
  //   glm::dvec3 shadowAtten = pLight->shadowAttenuation(r, p);
  //   glm::dvec3 attenuation = distanceAtten * shadowAtten;

  //   colorC += (diffuseTerm + specularTerm) * attenuation;
  // }
  // return glm::clamp(colorC, 0.0, 1.0);
  glm::dvec3 emissive = ke(i);        // Emissive term.
  glm::dvec3 ambient = ka(i);         // Ambient reflectance.
  glm::dvec3 diffuse = kd(i);         // Diffuse reflectance.
  glm::dvec3 specular = ks(i);        // Specular reflectance.
  double shininessVal = shininess(i); // Shininess exponent.

  // ADDED CODE: Get ambient light intensity from the scene.
  glm::dvec3 ambientLight = scene->ambient();

  // Compute the point of intersection.
  glm::dvec3 p = r.at(i);

  // Get the surface normal at the intersection (already interpolated).
  glm::dvec3 N = glm::normalize(i.getN());

  // Compute the view vector. (The ray r is cast from the eye to the point;
  // so the view direction is opposite to the ray direction.)
  glm::dvec3 V = glm::normalize(-r.getDirection());

  // Start with the emissive and ambient contributions.
  glm::dvec3 result = emissive + ambient * ambientLight;

  // ADDED CODE: Loop over all lights in the scene.
  const auto &lights = scene->getAllLights();
  for (const auto &light : lights)
  {
    // Compute the light direction (from the point to the light)
    glm::dvec3 L = glm::normalize(light->getDirection(p));

    // Compute the diffuse factor.
    // double diffFactor = glm::max(0.0, glm::dot(N, L));
    const Material& m = i.getMaterial();

    double diffFactor;
    if (m.Trans()) {
      // The material is transparent—process refraction.
       diffFactor = glm::abs(glm::dot(N, L));
  } else {
      diffFactor = glm::max(0.0, glm::dot(N, L));
  }
    //double diffFactor = glm::abs(glm::dot(N, L));
    // Compute the reflection direction. Here, the standard formula is used:
    // R = 2*(N dot L)*N - L.
    glm::dvec3 R = glm::normalize(2.0 * glm::dot(N, L) * N - L);

    // Compute the specular factor. Only add specular if the surface is lit.
    double specFactor = 0.0;
    if (diffFactor > 0.0)
      specFactor = pow(glm::max(0.0, glm::dot(R, V)), shininessVal);

    // Get the light color.
    glm::dvec3 lightColor = light->getColor();

    // Compute distance attenuation for this light.
    double distAtten = light->distanceAttenuation(p);

    // Compute shadow attenuation for this light.
    glm::dvec3 shadowAtten = light->shadowAttenuation(r, p);

    // Combine the diffuse and specular contributions.
    glm::dvec3 lightContribution = lightColor * distAtten * shadowAtten *
                                   (diffuse * diffFactor + specular * specFactor);

    // Add the light's contribution to the result.
    result += lightContribution;
  }

  // Optionally, clamp the final color components to the range [0, 1].
  result = glm::clamp(result, 0.0, 1.0);
  return result;
}

TextureMap::TextureMap(string filename) {
  data = readImage(filename.c_str(), width, height);
  if (data.empty()) {
    width = 0;
    height = 0;
    string error("Unable to load texture map '");
    error.append(filename);
    error.append("'.");
    throw TextureMapException(error);
  }
}

glm::dvec3 TextureMap::getMappedValue(const glm::dvec2 &coord) const {
  // YOUR CODE HERE
  //
  // In order to add texture mapping support to the
  // raytracer, you need to implement this function.
  // What this function should do is convert from
  // parametric space which is the unit square
  // [0, 1] x [0, 1] in 2-space to bitmap coordinates,
  // and use these to perform bilinear interpolation
  // of the values.
  double u = coord[0];
  double v = coord[1];
  
  // Map u,v to pixel space. Use (width-1) and (height-1) so that u==1 maps exactly to the last pixel.
  double uPos = u * (width - 1);
  double vPos = v * (height - 1);
  
  // Find the integer (floor) coordinates.
  int x = static_cast<int>(floor(uPos));
  int y = static_cast<int>(floor(vPos));
  
  // Compute the fractional part.
  double s = uPos - x;
  double t = vPos - y;
  
  // Clamp x and y so that x+1 and y+1 are in range.
  if(x < 0) x = 0;
  if(y < 0) y = 0;
  if(x >= width - 1) x = width - 2;
  if(y >= height - 1) y = height - 2;
  
  // Retrieve the four neighboring pixel values.
  glm::dvec3 c00 = getPixelAt(x, y);
  glm::dvec3 c10 = getPixelAt(x + 1, y);
  glm::dvec3 c01 = getPixelAt(x, y + 1);
  glm::dvec3 c11 = getPixelAt(x + 1, y + 1);
  
  // Bilinear interpolation:
  // First interpolate horizontally for the top and bottom rows.
  glm::dvec3 topInterp = (1.0 - s) * c00 + s * c10;
  glm::dvec3 bottomInterp = (1.0 - s) * c01 + s * c11;
  
  // Then interpolate vertically between the two results.
  glm::dvec3 finalColor = (1.0 - t) * topInterp + t * bottomInterp;
  
  return finalColor;
}

glm::dvec3 TextureMap::getPixelAt(int x, int y) const {
  // YOUR CODE HERE
  //
  // In order to add texture mapping support to the
  // raytracer, you need to implement this function.

  const uint8_t *pixel = data.data() + (x + y * getWidth()) * 3;
  double p = (double)pixel[0];
  double p1 = (double)pixel[1];
  double p2 = (double)pixel[2];
	return glm::dvec3(p / 255.0, p1 / 255.0, p2 / 255.0);

}

glm::dvec3 MaterialParameter::value(const isect &is) const {
  if (0 != _textureMap)
    return _textureMap->getMappedValue(is.getUVCoordinates());
  else
    return _value;
}

double MaterialParameter::intensityValue(const isect &is) const {
  if (0 != _textureMap) {
    glm::dvec3 value(_textureMap->getMappedValue(is.getUVCoordinates()));
    return (0.299 * value[0]) + (0.587 * value[1]) + (0.114 * value[2]);
  } else
    return (0.299 * _value[0]) + (0.587 * _value[1]) + (0.114 * _value[2]);
}
