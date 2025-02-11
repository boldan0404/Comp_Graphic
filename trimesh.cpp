#include "trimesh.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <float.h>
#include <string.h>
#include "../ui/TraceUI.h"
extern TraceUI *traceUI;
extern TraceUI *traceUI;

using namespace std;

Trimesh::~Trimesh()
{
  for (auto f : faces)
    delete f;
}

// must add vertices, normals, and materials IN ORDER
void Trimesh::addVertex(const glm::dvec3 &v) { vertices.emplace_back(v); }

void Trimesh::addNormal(const glm::dvec3 &n) { normals.emplace_back(n); }

void Trimesh::addColor(const glm::dvec3 &c) { vertColors.emplace_back(c); }

void Trimesh::addUV(const glm::dvec2 &uv) { uvCoords.emplace_back(uv); }

// Returns false if the vertices a,b,c don't all exist
bool Trimesh::addFace(int a, int b, int c)
{
  int vcnt = vertices.size();

  if (a >= vcnt || b >= vcnt || c >= vcnt)
    return false;

  TrimeshFace *newFace = new TrimeshFace(this, a, b, c);
  if (!newFace->degen)
    faces.push_back(newFace);
  else
    delete newFace;

  // Don't add faces to the scene's object list so we can cull by bounding
  // box
  return true;
}

// Check to make sure that if we have per-vertex materials or normals
// they are the right number.
const char *Trimesh::doubleCheck()
{
  if (!vertColors.empty() && vertColors.size() != vertices.size())
    return "Bad Trimesh: Wrong number of vertex colors.";
  if (!uvCoords.empty() && uvCoords.size() != vertices.size())
    return "Bad Trimesh: Wrong number of UV coordinates.";
  if (!normals.empty() && normals.size() != vertices.size())
    return "Bad Trimesh: Wrong number of normals.";

  return 0;
}

bool Trimesh::intersectLocal(ray &r, isect &i) const
{
  bool have_one = false;
  for (auto face : faces)
  {
    isect cur;
    if (face->intersectLocal(r, cur))
    {
      if (!have_one || (cur.getT() < i.getT()))
      {
        i = cur;
        have_one = true;
      }
    }
  }
  if (!have_one)
    i.setT(1000.0);
  return have_one;
}

bool TrimeshFace::intersect(ray &r, isect &i) const
{
  return intersectLocal(r, i);
}

// Intersect ray r with the triangle abc.  If it hits returns true,
// and put the parameter in t and the barycentric coordinates of the
// intersection in u (alpha) and v (beta).
bool TrimeshFace::intersectLocal(ray &r, isect &i) const
{
  // YOUR CODE HERE
  //
  // FIXME: Add ray-trimesh intersection

  /* To determine the color of an intersection, use the following rules:
     - If the parent mesh has non-empty `uvCoords`, barycentrically interpolate
       the UV coordinates of the three vertices of the face, then assign it to
       the intersection using i.setUVCoordinates().
     - Otherwise, if the parent mesh has non-empty `vertexColors`,
       barycentrically interpolate the colors from the three vertices of the
       face. Create a new material by copying the parent's material, set the
       diffuse color of this material to the interpolated color, and then
       assign this material to the intersection.
     - If neither is true, assign the parent's material to the intersection.
  */
  // Retrieve the triangle vertices from the parent mesh.
  const glm::dvec3 &V0 = parent->vertices[(*this)[0]];
  const glm::dvec3 &V1 = parent->vertices[(*this)[1]];
  const glm::dvec3 &V2 = parent->vertices[(*this)[2]];

  // Compute the edge vectors.
  glm::dvec3 E1 = V1 - V0;
  glm::dvec3 E2 = V2 - V0;

  // Compute the determinant using the cross product of the ray direction and E2.
  glm::dvec3 P = glm::cross(r.getDirection(), E2);
  double det = glm::dot(E1, P);

  // If the determinant is near zero, the ray is parallel to the triangle.
  if (det > -RAY_EPSILON && det < RAY_EPSILON)
    return false;

  double inv_det = 1.0 / det;
  glm::dvec3 T = r.getPosition() - V0;
  double alpha = glm::dot(T, P) * inv_det;

  // Check if the intersection point is outside the triangle.
  if (alpha < 0.0 || alpha > 1.0)
    return false;

  glm::dvec3 Q = glm::cross(T, E1);
  double beta = glm::dot(r.getDirection(), Q) * inv_det;

  if (beta < 0.0 || beta > 1.0 || (alpha + beta) > 1.0)
    return false;

  double t = glm::dot(E2, Q) * inv_det;
  if (t < RAY_EPSILON)
    return false;

  // Compute the third barycentric coordinate.
  double gamma = 1.0 - alpha - beta;

  // Set the intersection distance and object pointer.
  i.setT(t);
  i.setObject(this->parent);

  // ================================================================
  // ADDED CODE: Phong interpolation of normals.
  // If the parent mesh provides per-vertex normals, interpolate them
  // using barycentric coordinates and normalize the result.
  if (!parent->normals.empty())
  {
    glm::dvec3 normal0 = parent->normals[(*this)[0]];
    glm::dvec3 normal1 = parent->normals[(*this)[1]];
    glm::dvec3 normal2 = parent->normals[(*this)[2]];

    glm::dvec3 interpolatedNormal = glm::normalize(gamma * normal0 +
                                                   alpha * normal1 +
                                                   beta * normal2);

    i.setN(interpolatedNormal);
  }
  else
  {
    // Fallback to flat shading using the face normal.
    i.setN(glm::normalize(glm::cross(E1, E2)));
  }
  // ================================================================

  // ================================================================
  // ADDED CODE: Interpolation of per-vertex materials.
  // First, check for per-vertex UV coordinates. If they exist, use barycentric
  // interpolation to compute the UV coordinate at the intersection.
  if (!parent->uvCoords.empty())
  {
    glm::dvec2 uv0 = parent->uvCoords[(*this)[0]];
    glm::dvec2 uv1 = parent->uvCoords[(*this)[1]];
    glm::dvec2 uv2 = parent->uvCoords[(*this)[2]];
    glm::dvec2 interpolatedUV = gamma * uv0 + alpha * uv1 + beta * uv2;
    i.setUVCoordinates(interpolatedUV);
  }
  // Otherwise, if per-vertex colors (materials) exist, interpolate the diffuse
  // component without renormalization.
  else if (!parent->vertColors.empty())
  {
    glm::dvec3 color0 = parent->vertColors[(*this)[0]];
    glm::dvec3 color1 = parent->vertColors[(*this)[1]];
    glm::dvec3 color2 = parent->vertColors[(*this)[2]];
    glm::dvec3 interpolatedColor = gamma * color0 + alpha * color1 + beta * color2;
    Material mat = parent->getMaterial();
    mat.setDiffuse(interpolatedColor);
    i.setMaterial(mat);
  }
  // If no per-vertex material data exists, simply assign the parent's material.
  else
  {
    i.setMaterial(parent->getMaterial());
  }
  // ================================================================

  // Optionally, store the barycentric coordinates (if needed later).
  i.setBary(alpha, beta, gamma);

  return true;
  // i.setObject(this->parent);
  // return false;
}

// Once all the verts and faces are loaded, per vertex normals can be
// generated by averaging the normals of the neighboring faces.
void Trimesh::generateNormals()
{
  int cnt = vertices.size();
  normals.resize(cnt);
  std::vector<int> numFaces(cnt, 0);

  for (auto face : faces)
  {
    glm::dvec3 faceNormal = face->getNormal();

    for (int i = 0; i < 3; ++i)
    {
      normals[(*face)[i]] += faceNormal;
      ++numFaces[(*face)[i]];
    }
  }

  for (int i = 0; i < cnt; ++i)
  {
    if (numFaces[i])
      normals[i] /= numFaces[i];
  }

  vertNorms = true;
}