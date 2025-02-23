// kdTree.h
#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include "bbox.h"
#include "ray.h"

// A simple KD-tree template that accelerates ray–object intersection tests.
// T is expected to be a pointer type for scene objects (e.g., Geometry*)
// that provide a getBoundingBox() method and an intersect(ray&, isect&) method.
template <typename T>
class KdTree {
public:
    // Constructs a KD-tree from the given list of objects.
    // maxDepth: maximum recursion depth.
    // minObjects: minimum number of objects in a leaf node.
    KdTree(const std::vector<T>& objs, int maxDepth, int minObjects)
        : maxDepth(maxDepth), minObjects(minObjects)
    {
        // Make a local copy since we will partition it.
        std::vector<T> objects = objs;
        root = build(objects, 0);
    }
    
    ~KdTree() {
        destroy(root);
    }
    
    // Intersect a ray with the KD-tree.
    // Returns true if an intersection is found and updates 'closestIsect'
    // with the closest intersection (as determined by t).
    bool intersect(ray &r, isect &closestIsect) const {
        return intersectNode(root, r, closestIsect);
    }
    
private:
    // Internal KD-tree node.
    struct Node {
        BoundingBox bbox;       // Bounding box enclosing all objects in this node.
        Node* left;             // Left child.
        Node* right;            // Right child.
        std::vector<T> objects; // For leaf nodes, store the objects.
        bool isLeaf;            // True if this node is a leaf.
        
        Node() : left(nullptr), right(nullptr), isLeaf(false) { }
    };
    
    Node* root;
    int maxDepth;   // Maximum allowed depth.
    int minObjects; // Minimum number of objects per leaf.
    
    // Recursively builds the KD-tree using a simple median split.
    Node* build(std::vector<T>& objects, int depth) {
        Node* node = new Node();
        if (objects.empty())
            return node;
        
        // Compute the bounding box for this node.
        node->bbox = objects[0]->getBoundingBox();
        for (size_t i = 1; i < objects.size(); i++) {
            node->bbox.merge(objects[i]->getBoundingBox());
        }
        
        // Terminate recursion if node is small enough or max depth reached.
        if (depth >= maxDepth || objects.size() <= (size_t)minObjects) {
            node->isLeaf = true;
            node->objects = objects;
            return node;
        }
        
        // Choose the splitting axis: the one with the largest extent.
        glm::dvec3 extent = node->bbox.getMax() - node->bbox.getMin();
        int axis = 0;
        if (extent[1] > extent[0])
            axis = 1;
        if (extent[2] > extent[axis])
            axis = 2;
        
        // Sort objects by the center of their bounding box along the chosen axis.
        std::sort(objects.begin(), objects.end(), [axis](T a, T b) {
            glm::dvec3 centerA = (a->getBoundingBox().getMin() + a->getBoundingBox().getMax()) * 0.5;
            glm::dvec3 centerB = (b->getBoundingBox().getMin() + b->getBoundingBox().getMax()) * 0.5;
            return centerA[axis] < centerB[axis];
        });
        
        // Check the spread of centers along the chosen axis.
        glm::dvec3 centerFirst = (objects.front()->getBoundingBox().getMin() +
                                   objects.front()->getBoundingBox().getMax()) * 0.5;
        glm::dvec3 centerLast  = (objects.back()->getBoundingBox().getMin() +
                                   objects.back()->getBoundingBox().getMax()) * 0.5;
        if (fabs(centerLast[axis] - centerFirst[axis]) < 1e-6) {
            // The centers are nearly identical; do not split further.
            node->isLeaf = true;
            node->objects = objects;
            return node;
        }
        
        // Simple median split.
        size_t mid = objects.size() / 2;
        std::vector<T> leftObjs(objects.begin(), objects.begin() + mid);
        std::vector<T> rightObjs(objects.begin() + mid, objects.end());
        
        // If one side is empty, make this a leaf node.
        if (leftObjs.empty() || rightObjs.empty()) {
            node->isLeaf = true;
            node->objects = objects;
            return node;
        }
        
        // Recursively build child nodes.
        node->left = build(leftObjs, depth + 1);
        node->right = build(rightObjs, depth + 1);
        node->isLeaf = false;
        return node;
    }
    
    
    // Recursively tests the ray against the KD-tree nodes.
    bool intersectNode(Node* node, ray &r, isect &closestIsect) const {
        double tMin, tMax;
        if (!node || !node->bbox.intersect(r, tMin, tMax))
            return false;
        
        if (node->isLeaf) {
            bool hit = false;
            for (auto obj : node->objects) {
                isect tempIsect;
                if (obj->intersect(r, tempIsect)) {
                    if (tempIsect.getT() < closestIsect.getT()) {
                        closestIsect = tempIsect;
                        hit = true;
                    }
                }
            }
            return hit;
        }
        
        // Traverse both children.
        bool hitLeft = intersectNode(node->left, r, closestIsect);
        bool hitRight = intersectNode(node->right, r, closestIsect);
        return hitLeft || hitRight;
    }
    
    // Recursively deallocates nodes.
    void destroy(Node* node) {
        if (node) {
            destroy(node->left);
            destroy(node->right);
            delete node;
        }
    }
};

#endif // KDTREE_H
