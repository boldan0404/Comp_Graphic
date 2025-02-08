#include "cubeMap.h"
#include "../scene/material.h"
#include "../ui/TraceUI.h"
#include "ray.h"
extern TraceUI *traceUI;

glm::dvec3 CubeMap::getColor(ray r) const {
  // YOUR CODE HERE
  // FIXME: Implement Cube Map here
    // Normalize the ray direction.
    glm::dvec3 dir = glm::normalize(r.getDirection());

    // Compute absolute values for each component.
    double absX = fabs(dir.x);
    double absY = fabs(dir.y);
    double absZ = fabs(dir.z);

    int faceIndex;   // Which cube face to sample.
    double u, v;     // Texture coordinates in [0, 1].

    // Determine which face of the cube the ray direction points to.
    if (absX >= absY && absX >= absZ) {
        // X is the dominant component.
        if (dir.x > 0) {
            // Positive X face.
            faceIndex = 0; // tMap[0] is the +X face.
            u = 0.5 * (-dir.z / absX + 1.0);
            v = 0.5 * (-dir.y / absX + 1.0);
        } else {
            // Negative X face.
            faceIndex = 1; // tMap[1] is the -X face.
            u = 0.5 * (dir.z / absX + 1.0);
            v = 0.5 * (-dir.y / absX + 1.0);
        }
    } else if (absY >= absX && absY >= absZ) {
        // Y is the dominant component.
        if (dir.y > 0) {
            // Positive Y face.
            faceIndex = 2; // tMap[2] is the +Y face.
            u = 0.5 * (dir.x / absY + 1.0);
            v = 0.5 * (dir.z / absY + 1.0);
        } else {
            // Negative Y face.
            faceIndex = 3; // tMap[3] is the -Y face.
            u = 0.5 * (dir.x / absY + 1.0);
            v = 0.5 * (-dir.z / absY + 1.0);
        }
    } else {
        // Z is the dominant component.
        if (dir.z > 0) {
            // Positive Z face.
            faceIndex = 4; // tMap[4] is the +Z face.
            u = 0.5 * (dir.x / absZ + 1.0);
            v = 0.5 * (-dir.y / absZ + 1.0);
        } else {
            // Negative Z face.
            faceIndex = 5; // tMap[5] is the -Z face.
            u = 0.5 * (-dir.x / absZ + 1.0);
            v = 0.5 * (-dir.y / absZ + 1.0);
        }
    }

    // Now, if a texture map has been assigned for this face, use it.
    // Note: getMappedValue expects a glm::dvec2 in the [0, 1]×[0, 1] parametric space.
    if (tMap[faceIndex]) {
        return tMap[faceIndex]->getMappedValue(glm::dvec2(u, v));
    } else {
        // If no texture is assigned for this face, return black.
        return glm::dvec3(0.0, 0.0, 0.0);
    }
}

CubeMap::CubeMap() {}

CubeMap::~CubeMap() {}

void CubeMap::setNthMap(int n, TextureMap *m) {
  if (m != tMap[n].get())
    tMap[n].reset(m);
}
