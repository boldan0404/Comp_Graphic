#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

#define MAX_THREADS 32

// The main ray tracer.

#include "scene/cubeMap.h"
#include "scene/ray.h"
#include <glm/vec3.hpp>
#include <mutex>
#include <queue>
#include <thread>
#include <time.h>

class Scene;
class Pixel
{
public:
  Pixel(int i, int j, unsigned char *ptr) : ix(i), jy(j), value(ptr) {}

  int ix;
  int jy;
  unsigned char *value;
};

class RayTracer
{
public:
  RayTracer();
  ~RayTracer();

  glm::dvec3 tracePixel(int i, int j);
  glm::dvec3 traceRay(ray &r, const glm::dvec3 &thresh, int depth,
                      double &length);

  glm::dvec3 getPixel(int i, int j);
  void setPixel(int i, int j, glm::dvec3 color);
  void getBuffer(unsigned char *&buf, int &w, int &h);
  double aspectRatio();

  void traceImage(int w, int h);
  int aaImage();
  int RayTracer::adaptiveAntialiasImage(); //just added
  glm::dvec3 adaptiveSample(double x0, double y0, double x1, double y1, int depth);

  bool checkRender();
  void waitRender();

  void traceSetup(int w, int h);

  bool loadScene(const char *fn);
  bool sceneLoaded() { return scene != 0; }

  void setReady(bool ready) { m_bBufferReady = ready; }
  bool isReady() const { return m_bBufferReady; }
  bool refract(const glm::dvec3 &I, const glm::dvec3 &N, double eta, glm::dvec3 &T) const;
  const Scene &getScene() { return *scene; }

  bool stopTrace;

private:
  glm::dvec3 trace(double x, double y);

  std::unique_ptr<Scene> scene;
  std::vector<unsigned char> buffer;
  double thresh;
  int buffer_width, buffer_height;
  bool m_bBufferReady;

  int bufferSize;
  unsigned int threads;
  int block_size;
  double aaThresh;
  int samples;
};

#endif // __RAYTRACER_H__
