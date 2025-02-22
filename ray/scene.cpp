#include <cmath>

#include "../ui/TraceUI.h"
#include "kdTree.h"
#include "light.h"
#include "scene.h"
#include <glm/gtx/extended_min_max.hpp>
#include <glm/gtx/io.hpp>
#include <iostream>

using namespace std;

bool Geometry::intersect(ray &r, isect &i) const {
  double tmin, tmax;
  if (hasBoundingBoxCapability() && !(bounds.intersect(r, tmin, tmax)))
    return false;
  // Transform the ray into the object's local coordinate space
  glm::dvec3 pos = transform.globalToLocalCoords(r.getPosition());
  glm::dvec3 dir =
      transform.globalToLocalCoords(r.getPosition() + r.getDirection()) - pos;
  double length = glm::length(dir);
  dir = glm::normalize(dir);
  // Backup World pos/dir, and switch to local pos/dir
  glm::dvec3 Wpos = r.getPosition();
  glm::dvec3 Wdir = r.getDirection();
  r.setPosition(pos);
  r.setDirection(dir);
  bool rtrn = false;
  if (intersectLocal(r, i)) {
    // Transform the intersection point & normal returned back into
    // global space.
    i.setN(transform.localToGlobalCoordsNormal(i.getN()));
    i.setT(i.getT() / length);
    rtrn = true;
  }
  // Restore World pos/dir
  r.setPosition(Wpos);
  r.setDirection(Wdir);
  return rtrn;
}

bool Geometry::hasBoundingBoxCapability() const {
  // by default, primitives do not have to specify a bounding box. If this
  // method returns true for a primitive, then either the ComputeBoundingBox()
  // or the ComputeLocalBoundingBox() method must be implemented.

  // If no bounding box capability is supported for an object, that object will
  // be checked against every single ray drawn. This should be avoided whenever
  // possible, but this possibility exists so that new primitives will not have
  // to have bounding boxes implemented for them.
  return false;
}

void Geometry::ComputeBoundingBox() {
  // take the object's local bounding box, transform all 8 points on it,
  // and use those to find a new bounding box.

  BoundingBox localBounds = ComputeLocalBoundingBox();

  glm::dvec3 min = localBounds.getMin();
  glm::dvec3 max = localBounds.getMax();

  glm::dvec4 v, newMax, newMin;

  v = transform.localToGlobalCoords(glm::dvec4(min[0], min[1], min[2], 1));
  newMax = v;
  newMin = v;
  v = transform.localToGlobalCoords(glm::dvec4(max[0], min[1], min[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(min[0], max[1], min[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(max[0], max[1], min[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(min[0], min[1], max[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(max[0], min[1], max[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(min[0], max[1], max[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);
  v = transform.localToGlobalCoords(glm::dvec4(max[0], max[1], max[2], 1));
  newMax = glm::max(newMax, v);
  newMin = glm::min(newMin, v);

  bounds.setMax(glm::dvec3(newMax));
  bounds.setMin(glm::dvec3(newMin));
}

Scene::Scene() { ambientIntensity = glm::dvec3(0, 0, 0); }

Scene::~Scene() {
  for (auto &obj : objects)
    delete obj;
  for (auto &light : lights)
    delete light;
}

void Scene::buildKdTree(int treeDepth, int leafSize) {
  if (!objects.empty()) {
      kdtree = new KdTree<Geometry*>(objects, treeDepth, leafSize);
      if (kdtree != nullptr)
          std::cout << "KD-tree built with " << objects.size() << " objects." << std::endl;
      else
          std::cerr << "Error: KD-tree is null!" << std::endl;
  } else {
      kdtree = nullptr;
      std::cout << "No objects to build KD-tree." << std::endl;
  }
}



void Scene::add(Geometry *obj) {
  obj->ComputeBoundingBox();
  sceneBounds.merge(obj->getBoundingBox());
  objects.emplace_back(obj);
}

void Scene::add(Light *light) { lights.emplace_back(light); }


// Get any intersection with an object.  Return information about the
// intersection through the reference parameter.
bool Scene::intersect(ray &r, isect &i) const {
  if (kdtree) {
      // Initialize 'i' with a very large t value.
      i.setT(1.0e308);
      if (kdtree->intersect(r, i))
          return true;
      return false;
  }
  // Fallback: linear search over objects
  bool hit = false;
  for (const auto &obj : objects) {
      isect temp;
      if (obj->intersect(r, temp)) {
          if (!hit || (temp.getT() < i.getT())) {
              i = temp;
              hit = true;
          }
      }
  }
  return hit;
}

TextureMap *Scene::getTexture(string name) {
  auto itr = textureCache.find(name);
  if (itr == textureCache.end()) {
    textureCache[name].reset(new TextureMap(name));
    return textureCache[name].get();
  }
  return itr->second.get();
}
