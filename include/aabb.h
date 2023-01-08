#ifndef AABB_H_
#define AABB_H_
#include "core.h"
#include "ray.h"

class AABB {
	[[nodiscard]] virtual int maxExtent() const = 0;
	[[nodiscard]] virtual float SurfaceArea() const = 0;
};

class Bounds2 : public AABB {
public:
	Vec2f pMin, pMax;// two points to specify the bounding box
	Bounds2();
	Bounds2(const Vec2f &p1, const Vec2f &p2);
	Bounds2(const float &uMin, const float &vMin, const float &uMax, const float &vMax);
	[[nodiscard]] int maxExtent() const override;
	[[nodiscard]] Vec2f Diagonal() const;
	[[nodiscard]] float SurfaceArea() const override;
	[[nodiscard]] Vec2f Centroid() const;
	[[nodiscard]] Bounds2 Intersect(const Bounds2 &b) const;
	[[nodiscard]] Vec2f Offset(const Vec2f &p) const;
	[[nodiscard]] static bool Overlaps(const Bounds2 &b1, const Bounds2 &b2);
	[[nodiscard]] static bool Inside(const Vec2f &p, const Bounds2 &b);
	[[nodiscard]] bool Inside(const Vec2f &p);
};

class Bounds3 : public AABB {
public:
	Vec3f pMin, pMax;// two points to specify the bounding box
	Bounds3();
	explicit Bounds3(const Vec3f &p);
	Bounds3(const Vec3f &p1, const Vec3f &p2);
	Bounds3(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3);
	[[nodiscard]] Vec3f Diagonal() const;
	[[nodiscard]] int maxExtent() const override;
	[[nodiscard]] float SurfaceArea() const override;
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

Bounds2 Union(const Bounds2 &b1, const Bounds2 &b2);
Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2);
Bounds3 Union(const Bounds3 &b, const Vec3f &p);
#endif// AABB_H_