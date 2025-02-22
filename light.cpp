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
  glm::dvec3 atten = r.getAtten();
  glm::dvec3 lightDir = glm::normalize(getDirection(p));

  glm::dvec3 orgin = p + RAY_EPSILON * lightDir;
  // Create a shadow ray starting at p (offset by RAY_EPSILON).
  ray shadowRay(orgin, lightDir, atten, ray::SHADOW);
  isect i;

  // Traverse the ray until no more intersections are found, or an intersection occurs beyond the light.
  if (scene->intersect(shadowRay, i))
  {
    const Material &m = i.getMaterial();

    if (!m.Trans())
    {
      // if its opaque
      return glm::dvec3(0.0, 0.0, 0.0);
    }
    else
    {
      // if its transparent
      // doucument the ray enter point
      shadowRay.setPosition(shadowRay.at(i) + RAY_EPSILON * lightDir);

      isect exitI;
      if (scene->intersect(shadowRay, exitI))
      {
        double exitT = exitI.getT();
        atten *= glm::pow(m.kt(i), glm::dvec3(exitT));
      }
      else
      {
        double d = i.getT();
        atten *= glm::pow(m.kt(i), glm::dvec3(d));
      }

      // shadowRay = ray(shadowRay.at(exitI) + RAY_EPSILON * lightDir, lightDir, atten, ray::SHADOW);
    }
  }
  return atten;

  //-----Jerry code-----
  // glm::dvec3 lightDir = glm::normalize(r.getDirection(p)); // This is defined below as -orientation.

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
  double denominator = constantTerm + linearTerm * d + quadraticTerm * d * d;
  double attenuation = 1.0 / denominator;
  if (attenuation > 1.0)
  {
    attenuation = 1.0;
  }
  return attenuation;
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
  // // --------Ann  code-------
  glm::dvec3 atten = r.getAtten();
  glm::dvec3 lightDir = glm::normalize(getDirection(p));
  glm::dvec3 orgin = p + RAY_EPSILON * lightDir;
  // Create a shadow ray starting at p (offset by RAY_EPSILON).
  ray shadowRay(orgin, lightDir, atten, ray::SHADOW);
  isect i;
  // For a point light, we need to consider only intersections between p and the light.
  double maxDistance = glm::distance(p, position);

  // Traverse the ray until no more intersections are found, or an intersection occurs beyond the light.
  if (scene->intersect(shadowRay, i))
  {
    // If the hit is farther than the light, then no (further) occluders matter.
    if (i.getT() > maxDistance)
      return atten;

    const Material &m = i.getMaterial();

    if (!m.Trans())
    {
      // if its opaque
      return glm::dvec3(0.0, 0.0, 0.0);
    }
    else
    {
      // if its transparent
      // doucument the ray enter point
      shadowRay.setPosition(shadowRay.at(i) + RAY_EPSILON * lightDir);

      isect exitI;
      if (scene->intersect(shadowRay, exitI))
      {
        double exitT = exitI.getT();
        atten *= glm::pow(m.kt(i), glm::dvec3(exitT));
      }
      else
      {
        double d = i.getT();
        atten *= glm::pow(m.kt(i), glm::dvec3(d));
      }
      // shadowRay = ray(shadowRay.at(exitI) + RAY_EPSILON * lightDir, lightDir, atten, ray::SHADOW);
    }
  }
  return atten;

  // //------Jerry code-----
  // glm::dvec3 lightDir = glm::normalize(r.getDirection(p));
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
