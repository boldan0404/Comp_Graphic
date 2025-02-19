#include "cubeMap.h"
#include "../scene/material.h"
#include "../ui/TraceUI.h"
#include "ray.h"
extern TraceUI *traceUI;

glm::dvec3 CubeMap::getColor(ray r) const
{
    // YOUR CODE HERE
    // FIXME: Implement Cube Map here
    glm::dvec3 dir = r.getDirection();
    int face;
    double u, v;

    double absX = std::abs(dir.x);
    double absY = std::abs(dir.y);
    double absZ = std::abs(dir.z);

    // Determine which face of the cube map to sample.
    if ((absX > absY) && (absX > absZ))
    {
        if (dir.x > 0)
        {
            face = 0; // +X
            u = (dir.z / absX + 1.0) * 0.5;
            v = (dir.y / absX + 1.0) * 0.5;
        }
        else
        {
            face = 1; // -X
            u = (-dir.z / absX + 1.0) * 0.5;
            v = (dir.y / absX + 1.0) * 0.5;
        }
    }
    else if ((absY > absX) && (absY > absZ))
    {
        if (dir.y > 0)
        {
            face = 2; // +Y
            u = (dir.x / absY + 1.0) * 0.5;
            v = (dir.z / absY + 1.0) * 0.5;
        }
        else
        {
            face = 3; // -Y
            u = (dir.x / absY + 1.0) * 0.5;
            v = (-dir.z / absY + 1.0) * 0.5;
        }
    }
    else if ((absZ > absY) && (absZ > absX))
    {
        if (dir.z > 0)
        {
            face = 5; // +Z
            u = (-dir.x / absZ + 1.0) * 0.5;
            v = (dir.y / absZ + 1.0) * 0.5;
        }
        else
        {
            face = 4; // -Z
            u = (dir.x / absZ + 1.0) * 0.5;
            v = (dir.y / absZ + 1.0) * 0.5;
        }
    }

    // Sample the cube map.
    if (tMap[face])
    {
        return tMap[face]->getMappedValue(glm::dvec2(u, v));
    }
    else
    {
        return glm::dvec3(0, 0, 0);
    }
}

CubeMap::CubeMap() {}

CubeMap::~CubeMap() {}

void CubeMap::setNthMap(int n, TextureMap *m)
{
    if (m != tMap[n].get())
        tMap[n].reset(m);
}