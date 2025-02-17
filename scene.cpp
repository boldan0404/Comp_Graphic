#include <cmath>

#include "../ui/TraceUI.h"
#include "kdTree.h"
#include "light.h"
#include "scene.h"
#include <glm/gtx/extended_min_max.hpp>
#include <glm/gtx/io.hpp>
#include <iostream>

using namespace std;

Scene::BVHNode* Scene::buildBVH(std::vector<Geometry*>& objs, int depth) {
  BVHNode* node = new BVHNode();

  // Compute the bounding box for all objects in objs.
  BoundingBox nodeBounds;
  bool first = true;
  for (auto obj : objs) {
    if (first) {
      nodeBounds = obj->getBoundingBox();
      first = false;
    } else {
      nodeBounds.merge(obj->getBoundingBox());
    }
  }
  node->bounds = nodeBounds;

  // If few objects remain, make this a leaf node.
  if (objs.size() <= 4) {
    node->primitives = objs;
    return node;
}

  // Choose splitting axis: axis with largest extent.
  glm::dvec3 extent = nodeBounds.getMax() - nodeBounds.getMin();
  int axis = 0;
  if (extent[1] > extent[0])
      axis = 1;
  if (extent[2] > extent[axis])
      axis = 2;

  // Sort objects by the center of their bounding boxes along the chosen axis.
  std::sort(objs.begin(), objs.end(), [axis](Geometry* a, Geometry* b) {
      double aCenter = (a->getBoundingBox().getMin()[axis] + a->getBoundingBox().getMax()[axis]) * 0.5;
      double bCenter = (b->getBoundingBox().getMin()[axis] + b->getBoundingBox().getMax()[axis]) * 0.5;
      return aCenter < bCenter;
  });

  // Split objects into two halves.
  size_t mid = objs.size() / 2;
  std::vector<Geometry*> leftObjs(objs.begin(), objs.begin() + mid);
  std::vector<Geometry*> rightObjs(objs.begin() + mid, objs.end());

  node->left.reset(buildBVH(leftObjs, depth + 1));
  node->right.reset(buildBVH(rightObjs, depth + 1));

  return node;
}

void Scene::buildBVH() {
  // Build the BVH using all objects in the scene.
  std::vector<Geometry*> objs = objects;  // Make a copy so we can sort.
  bvhRoot.reset(buildBVH(objs, 0));
}

// ---------------- BVH Intersection ----------------

bool Scene::intersectBVH(BVHNode* node, ray &r, isect &hitRecord) const {
  double tmin, tmax;
  // If the ray does not hit this node's bounding box, then no intersections occur here.
  if (!node->bounds.intersect(r, tmin, tmax))
    return false;

  // If this is a leaf node (no children), test all primitives.
  if (!node->left && !node->right) {
    bool hitSomething = false;
    double closestT = std::numeric_limits<double>::max();
    isect tempRecord;
    for (auto obj : node->primitives) {
      if (obj->intersect(r, tempRecord)) {
        if (tempRecord.getT() < closestT) {
          closestT = tempRecord.getT();
          hitRecord = tempRecord;
          hitSomething = true;
        }
      }
    }
    return hitSomething;
  }

  // Otherwise, this is an internal node.
  bool hitLeft = false, hitRight = false;
  isect leftRecord, rightRecord;
  if (node->left)
    hitLeft = intersectBVH(node->left.get(), r, leftRecord);
  if (node->right)
    hitRight = intersectBVH(node->right.get(), r, rightRecord);

  // If both children are hit, choose the closer intersection.
  if (hitLeft && hitRight) {
    if (leftRecord.getT() < rightRecord.getT())
      hitRecord = leftRecord;
    else
      hitRecord = rightRecord;
    return true;
  } else if (hitLeft) {
    hitRecord = leftRecord;
    return true;
  } else if (hitRight) {
    hitRecord = rightRecord;
    return true;
  }

  return false;
}


// ---------------- Modified Scene::intersect ----------------

bool Scene::intersect(ray &r, isect &i) const {
  if (bvhRoot) {
    bool hit = intersectBVH(bvhRoot.get(), r, i);
    if (!hit)
      i.setT(1000.0);
    // if (TraceUI::m_debug) {
    //   addToIntersectCache(std::make_pair(new ray(r), new isect(i)));
    // }
    return hit;
  }
  // Fallback: linear search over all objects.
  bool hitAnything = false;
  for (const auto &obj : objects) {
    isect cur;
    if (obj->intersect(r, cur)) {
      if (!hitAnything || (cur.getT() < i.getT())) {
        i = cur;
        hitAnything = true;
      }
    }
  }
  if (!hitAnything)
    i.setT(1000.0);
  // if (TraceUI::m_debug) {
  //   addToIntersectCache(std::make_pair(new ray(r), new isect(i)));
  // }
  return hitAnything;
}

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

void Scene::add(Geometry *obj) {
  obj->ComputeBoundingBox();
  sceneBounds.merge(obj->getBoundingBox());
  objects.emplace_back(obj);
}

void Scene::add(Light *light) { lights.emplace_back(light); }


// Get any intersection with an object.  Return information about the
// intersection through the reference parameter.
// bool Scene::intersect(ray &r, isect &i) const {
//   double tmin = 0.0;
//   double tmax = 0.0;
//   bool have_one = false;
//   for (const auto &obj : objects) {
//     isect cur;
//     if (obj->intersect(r, cur)) {
//       if (!have_one || (cur.getT() < i.getT())) {
//         i = cur;
//         have_one = true;
//       }
//     }
//   }
//   if (!have_one)
//     i.setT(1000.0);
//   // if debugging,
//   if (TraceUI::m_debug) {
//     addToIntersectCache(std::make_pair(new ray(r), new isect(i)));
//   }
//   return have_one;
// }

TextureMap *Scene::getTexture(string name) {
  auto itr = textureCache.find(name);
  if (itr == textureCache.end()) {
    textureCache[name].reset(new TextureMap(name));
    return textureCache[name].get();
  }
  return itr->second.get();
}
