#include <cmath>
#include <iostream>

#include "light.h"
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>

using namespace std;

double DirectionalLight::distanceAttenuation(const glm::dvec3 &) const
{
  // distance to light is infinite, so f(di) goes to 0.  Return 1.
  return 1.0;
}

glm::dvec3 DirectionalLight::shadowAttenuation(const ray &r,
                                               const glm::dvec3 &p) const
{
  // YOUR CODE HERE:
  // You should implement shadow-handling code here.
  //-----Ann code-----
  // create a direction from intersection point to light source
  glm::dvec3 lightDir = getDirection(p);

  // create a shadow ray from intersection point to light source(where do i get the intersection point?)
  ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1.0, 1.0, 1.0), ray::SHADOW);
  isect i;

  glm::dvec3 I_before(1.0, 1.0, 1.0);
  // check if the shadow ray intersects with any object
  while (scene->intersect(shadowRay, i))
  {
    // check material property of the object
    const Material &m = i.getMaterial();

    if (m.kt(i) == glm::dvec3(0.0, 0.0, 0.0))
    {
      // if its opaque
      return glm::dvec3(0.0, 0.0, 0.0);
    }
    else
    {
      // for translucent object
      double enrtyT = i.getT();
      glm::dvec3 entryPoint = shadowRay.at(i);

      ray exitRay(entryPoint + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
      isect exitIntersection;

      double exitT = exitIntersection.getT();
      double d = exitT - enrtyT; // Distance light travels through the object

      glm::dvec3 transparency = m.kt(i);
      I_before *= glm::pow(transparency, glm::dvec3(d, d, d));

      shadowRay = ray(exitRay.at(exitIntersection) + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
    }
  }
  // if it does not, return the color of the light source
  return I_before;

  //-----Jerry code-----
  // glm::dvec3 lightDir = getDirection(p); // This is defined below as -orientation.

  // // Construct a shadow ray starting just above p.
  // ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);

  // isect i;
  // // Start with full light coming in.
  // glm::dvec3 transmission(1.0, 1.0, 1.0);

  // // For directional lights the shadow ray goes to infinity.
  // // We repeatedly check for intersections along the ray; if an intersected object is transparent,
  // // we multiply the shadow attenuation by its transparency (kt). Otherwise, if an opaque object is hit,
  // // the final transmission is zero.
  // while (scene->intersect(shadowRay, i))
  // {
  //   transmission *= i.getMaterial().kt(i);

  //   // If the material is completely opaque (or transmission has dropped very low) then we can quit.
  //   if (transmission[0] <= 0.0 && transmission[1] <= 0.0 && transmission[2] <= 0.0)
  //     return glm::dvec3(0.0, 0.0, 0.0);

  //   // Advance the shadow ray origin to just past the current intersection.
  //   shadowRay = ray(shadowRay.at(i) + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
  // }
  // return transmission;
}

glm::dvec3 DirectionalLight::getColor() const { return color; }

glm::dvec3 DirectionalLight::getDirection(const glm::dvec3 &) const
{
  return -orientation;
}

double PointLight::distanceAttenuation(const glm::dvec3 &P) const
{
  // YOUR CODE HERE
  double d = glm::distance(P, position);

  double constantTerm = 1.0;
  double linearTerm = 0.001;
  double quadraticTerm = 0.0001;
  double attenuation = 1.0 / (constantTerm + linearTerm * d + quadraticTerm * d * d);

  if (attenuation > 1.0)
    attenuation = 1.0;

  return std::min(attenuation, 1.0);
}

glm::dvec3 PointLight::getColor() const { return color; }

glm::dvec3 PointLight::getDirection(const glm::dvec3 &P) const
{
  return glm::normalize(position - P);
}

glm::dvec3 PointLight::shadowAttenuation(const ray &r,
                                         const glm::dvec3 &p) const
{
  // YOUR CODE HERE:
  // You should implement shadow-handling code here.
  // --------Ann  code-------
  // Compute the (normalized) direction from point p to the light.
  glm::dvec3 lightDir = getDirection(p);

  // Create a shadow ray starting at p (offset by RAY_EPSILON).
  ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
  isect i;

  // Begin with full light coming in.
  glm::dvec3 I_before(1.0, 1.0, 1.0);

  // For a point light, we need to consider only intersections between p and the light.
  double maxDistance = glm::distance(p, position);

  // Traverse the ray until no more intersections are found, or an intersection occurs beyond the light.
  while (scene->intersect(shadowRay, i))
  {
    // If the hit is farther than the light, then no (further) occluders matter.
    glm::dvec3 entryPoint = shadowRay.at(i);
    double entryT = i.getT();
    if (entryT > maxDistance)
      break;

    // Move the ray origin slightly inside to find the exit point
    ray exitRay(entryPoint + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
    isect exitIntersection;

    double exitT = exitIntersection.getT();
    double d = exitT - entryT; // Distance light travels through the object

    // calculate material transparency
    glm::dvec3 transparency = i.getMaterial().kt(i);
    I_before *= glm::pow(transparency, glm::dvec3(d, d, d));

    // If the cumulative transparency is zero, we can return immediately.
    if (I_before[0] <= 0.0 && I_before[1] <= 0.0 && I_before[2] <= 0.0)
      return glm::dvec3(0.0, 0.0, 0.0);

    // Move the ray origin to just past the exit point to check for more objects
    shadowRay = ray(exitRay.at(exitIntersection) + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
  }
  return I_before;

  //------Jerry code-----
  // glm::dvec3 lightDir = getDirection(p);
  // // Create a shadow ray starting at p (offset by RAY_EPSILON).
  // ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
  // isect i;

  // // Begin with full light coming in.
  // glm::dvec3 transmission(1.0, 1.0, 1.0);

  // // For a point light, we need to consider only intersections between p and the light.
  // double maxDistance = glm::distance(p, position);

  // // Traverse the ray until no more intersections are found, or an intersection occurs beyond the light.
  // while (scene->intersect(shadowRay, i))
  // {
  //   // If the hit is farther than the light, then no (further) occluders matter.
  //   if (i.getT() > maxDistance)
  //     break;

  //   // Multiply the current transmission by the intersected material's transparency.
  //   transmission *= i.getMaterial().kt(i);

  //   // If the cumulative transparency is zero, we can return immediately.
  //   if (transmission[0] <= 0.0 && transmission[1] <= 0.0 && transmission[2] <= 0.0)
  //     return glm::dvec3(0.0, 0.0, 0.0);

  //   // Advance the shadow ray past the current intersection.
  //   shadowRay = ray(shadowRay.at(i) + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1, 1, 1), ray::SHADOW);
  // }
  // return transmission;
}

#define VERBOSE 0
