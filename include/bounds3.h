#ifndef BOUNDS3_H_
#define BOUNDS3_H_
#include "core.h"
#include "ray.h"

class Bounds3 {
public:
	Vec3f pMin, pMax;// two points to specify the bounding box
	Bounds3();
	explicit Bounds3(const Vec3f &p);
	Bounds3(const Vec3f &p1, const Vec3f &p2);
	Bounds3(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3);
	[[nodiscard]] Vec3f Diagonal() const;
	[[nodiscard]] int maxExtent() const;
	[[nodiscard]] float SurfaceArea() const;
	[[nodiscard]] Vec3f Centroid() const;
	[[nodiscard]] Bounds3 Intersect(const Bounds3 &b) const;
	[[nodiscard]] Vec3f Offset(const Vec3f &p) const;
	[[nodiscard]] static bool Overlaps(const Bounds3 &b1, const Bounds3 &b2);
	[[nodiscard]] static bool Inside(const Vec3f &p, const Bounds3 &b);
	const Vec3f &operator[](int i) const;
	[[nodiscard]] bool IntersectP(const Ray &ray) const;
	bool IntersectP(const Ray &ray, float &tMin, float &tMax) const;
	[[nodiscard]] std::pair<bool, float> IntersectWithTime(const Ray &ray) const;
};

Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2);
Bounds3 Union(const Bounds3 &b, const Vec3f &p);
#endif// BOUNDS3_H_