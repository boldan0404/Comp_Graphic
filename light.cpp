#include <cmath>
#include <iostream>

#include "light.h"
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>

using namespace std;

double DirectionalLight::distanceAttenuation(const glm::dvec3 &) const {
  // distance to light is infinite, so f(di) goes to 0.  Return 1.
  return 1.0;
}

// Revised version for directional light shadow attenuation (recursive version)
glm::dvec3 DirectionalLight::shadowAttenuation(const ray &r, const glm::dvec3 &p) const {
  // Compute the constant light direction (for a directional light, getDirection returns -orientation)
  glm::dvec3 lightDir = glm::normalize(getDirection(p));

  // Construct the initial shadow ray starting slightly above p
  ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, r.getAtten(), ray::SHADOW);

  // Record the first intersection along the shadow ray
  isect firstHit;
  if (scene->intersect(shadowRay, firstHit)) {
    // If an intersection is found, check if the material is transparent
    if (firstHit.getMaterial().Trans()) {
      // Get the parametric distance to the first hit
      double t_first = firstHit.getT();

      // Save the current ray's direction and attenuation for clarity
      glm::dvec3 currentDir = shadowRay.getDirection();
      glm::dvec3 currentAtten = shadowRay.getAtten();

      // Create a new shadow ray starting just beyond the first intersection
      ray nextRay(shadowRay.at(t_first + RAY_EPSILON), glm::normalize(currentDir), currentAtten, ray::SHADOW);

      // Record the next intersection along this new ray
      isect secondHit;
      if (scene->intersect(nextRay, secondHit)) {
        double t_second = secondHit.getT();
        // Compute the distance between the first hit point and the second hit point
        double segmentLength = glm::distance(shadowRay.at(firstHit), nextRay.at(secondHit));
        // Recursively compute the attenuation from just past the second hit
        glm::dvec3 remainingAtten = shadowAttenuation(nextRay, nextRay.at(t_second + RAY_EPSILON));
        // Multiply the local attenuation (based on the first hit) by the remaining attenuation
        return glm::pow(firstHit.getMaterial().kt(firstHit), glm::dvec3(segmentLength)) * remainingAtten;
      }
      else {
        // If no further intersection is found, compute attenuation from p to the first hit
        double segmentLength = glm::distance(p, shadowRay.at(firstHit));
        return glm::pow(firstHit.getMaterial().kt(firstHit), glm::dvec3(segmentLength)) * getColor();
      }
    }
    // If the material is opaque, return zero light
    return glm::dvec3(0.0, 0.0, 0.0);
  }
  else {
    // If no intersection is found along the shadow ray, return full light
    return getColor();
  }
}


glm::dvec3 DirectionalLight::getColor() const { return color; }

glm::dvec3 DirectionalLight::getDirection(const glm::dvec3 &) const {
  return -orientation;
}

double PointLight::distanceAttenuation(const glm::dvec3 &P) const {
  // YOUR CODE HERE

  // You'll need to modify this method to attenuate the intensity
  // of the light based on the distance between the source and the
  // point P.  For now, we assume no attenuation and just return 1.0

    double d = glm::distance(P, position);
    // Attenuation coefficients (tweak these as needed)
    // double constantTerm  = 1.0;
    // double linearTerm    = 0.001;
    // double quadraticTerm = 0.0001;
    double denominator = constantTerm + linearTerm * d + quadraticTerm * d * d;
    double attenuation = 1.0 / denominator;
    if (attenuation > 1.0)
    {
       attenuation = 1.0;
    }
       return attenuation;
}

glm::dvec3 PointLight::getColor() const { return color; }

glm::dvec3 PointLight::getDirection(const glm::dvec3 &P) const {
  return glm::normalize(position - P);
}

glm::dvec3 PointLight::shadowAttenuation(const ray &r, const glm::dvec3 &p) const {
  // Compute the light direction from p to the point light.
  glm::dvec3 lightDir = glm::normalize(getDirection(p));  // getDirection(p) returns (position - p) normalized.
  
  // Compute the maximum distance from the shading point to the light.
  double maxDistance = glm::distance(p, position);
  
  // Construct the initial shadow ray starting just off p.
  ray shadowRay(p + RAY_EPSILON * lightDir, lightDir, glm::dvec3(1,1,1), ray::SHADOW);
  
  // Record the first intersection along the shadow ray.
  isect firstHit;
  if (scene->intersect(shadowRay, firstHit) && firstHit.getT() <= maxDistance) {
    // If an intersection is found and it occurs before reaching the light,
    // check if the intersected material is transparent.
    if (firstHit.getMaterial().Trans()) {
      // Get the distance along the shadow ray to the first hit.
      double t_first = firstHit.getT();
      
      // Create a new shadow ray starting just beyond the first hit.
      ray nextRay(shadowRay.at(t_first + RAY_EPSILON), glm::normalize(shadowRay.getDirection()), shadowRay.getAtten(), ray::SHADOW);
      
      // Record the next intersection along the new shadow ray.
      isect secondHit;
      if (scene->intersect(nextRay, secondHit)) {
        double t_second = secondHit.getT();
        // Compute the distance between the first hit point and this next hit point.
        double segmentLength = glm::distance(shadowRay.at(firstHit), nextRay.at(secondHit));
        // Recursively compute the attenuation from just beyond the second hit.
        glm::dvec3 remainingAtten = shadowAttenuation(nextRay, nextRay.at(t_second + RAY_EPSILON));
        // Return the product of the local attenuation and the remaining attenuation.
        return glm::pow(firstHit.getMaterial().kt(firstHit), glm::dvec3(segmentLength)) * remainingAtten;
      }
      else {
        // If no further intersection is found, compute attenuation from p to the first hit.
        double segmentLength = glm::distance(p, shadowRay.at(firstHit));
        return glm::pow(firstHit.getMaterial().kt(firstHit), glm::dvec3(segmentLength)) * getColor();
      }
    }
    // If the intersected material is opaque, return zero light.
    return glm::dvec3(0.0, 0.0, 0.0);
  }
  else {
    // If no intersection is found or the hit is beyond the light, return full light.
    return getColor();
  }
}


#define VERBOSE 0

