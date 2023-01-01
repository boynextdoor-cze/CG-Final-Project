#ifndef BOUNDS3_H_
#define BOUNDS3_H_
#include "core.h"
#include "ray.h"

class Bounds3
{
public:
  Vec3f pMin, pMax; // two points to specify the bounding box
  Bounds3();
  Bounds3(const Vec3f p);
  Bounds3(const Vec3f p1, const Vec3f p2);
  Bounds3(const Vec3f p1, const Vec3f p2, const Vec3f p3);
  Vec3f Diagonal() const;
  int maxExtent() const;
  float SurfaceArea() const;
  Vec3f Centroid() const;
  Bounds3 Intersect(const Bounds3 &b) const;
  Vec3f Offset(const Vec3f &p) const;
  bool Overlaps(const Bounds3 &b1, const Bounds3 &b2) const;
  bool Inside(const Vec3f &p, const Bounds3 &b) const;
  const Vec3f &operator[](int i) const;
  bool IntersectP(const Ray &ray) const;
	bool IntersectP(const Ray &ray, float &tMin, float &tMax) const;
  std::pair<bool, float> IntersectWithTime(const Ray &ray) const;
};

Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2);
Bounds3 Union(const Bounds3 &b, const Vec3f &p);
#endif //BOUND3_H_ 