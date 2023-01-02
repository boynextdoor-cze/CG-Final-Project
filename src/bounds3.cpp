#include "bounds3.h"

const float float_min = std::numeric_limits<float>::lowest();
const float float_max = std::numeric_limits<float>::max();

Bounds3::Bounds3() {
	pMax = Vec3f(float_min, float_min, float_min);
	pMin = Vec3f(float_max, float_max, float_max);
}

Bounds3::Bounds3(const Vec3f &p) : pMin(p), pMax(p) {}

Bounds3::Bounds3(const Vec3f &p1, const Vec3f &p2) {
	pMin = p1.cwiseMin(p2);
	pMax = p1.cwiseMax(p2);
}

Bounds3::Bounds3(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
	pMin = p1.cwiseMin(p2.cwiseMin(p3));
	pMax = p1.cwiseMax(p2.cwiseMax(p3));
}

Vec3f Bounds3::Diagonal() const { return pMax - pMin; }

int Bounds3::maxExtent() const {
	Vec3f d = Diagonal();
	if (d.x() > d.y() && d.x() > d.z())
		return 0;
	else if (d.y() > d.z())
		return 1;
	else
		return 2;
}

float Bounds3::SurfaceArea() const {
	Vec3f d = Diagonal();
	return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
}

Vec3f Bounds3::Centroid() const { return 0.5 * pMin + 0.5 * pMax; }

Bounds3 Bounds3::Intersect(const Bounds3 &b) const {
	return Bounds3{pMin.cwiseMax(b.pMin), pMax.cwiseMin(b.pMax)};
}

Vec3f Bounds3::Offset(const Vec3f &p) const {
	Vec3f o = p - pMin;
	if (pMax.x() > pMin.x())
		o.x() /= pMax.x() - pMin.x();
	if (pMax.y() > pMin.y())
		o.y() /= pMax.y() - pMin.y();
	if (pMax.z() > pMin.z())
		o.z() /= pMax.z() - pMin.z();
	return o;
}

bool Bounds3::Overlaps(const Bounds3 &b1, const Bounds3 &b2) {
	bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
	bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
	bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
	return (x && y && z);
}

bool Bounds3::Inside(const Vec3f &p, const Bounds3 &b) {
	return (p.x() >= b.pMin.x() - EPS && p.x() <= b.pMax.x() + EPS &&
	        p.y() >= b.pMin.y() - EPS && p.y() <= b.pMax.y() + EPS &&
	        p.z() >= b.pMin.z() - EPS && p.z() <= b.pMax.z() + EPS);
}
const Vec3f &Bounds3::operator[](int i) const { return (i == 0) ? pMin : pMax; }

bool Bounds3::IntersectP(const Ray &ray) const {
	// invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because
	// Multiply is faster that Division
	float t_l = float_min, t_r = float_max;
	Vec3f l = (pMin - ray.origin).cwiseProduct(ray.inv_direction);
	Vec3f r = (pMax - ray.origin).cwiseProduct(ray.inv_direction);
	t_l = l.cwiseMin(r).maxCoeff();
	t_r = l.cwiseMax(r).minCoeff();
	if (t_l > t_r + EPS)
		return false;
	if (t_r < -EPS)
		return false;
	return true;
}

bool Bounds3::IntersectP(const Ray &ray, float &tMin, float &tMax) const {
	// invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because
	// Multiply is faster that Division
	float t_l = float_min, t_r = float_max;
	Vec3f l = (pMin - ray.origin).cwiseProduct(ray.inv_direction);
	Vec3f r = (pMax - ray.origin).cwiseProduct(ray.inv_direction);
	t_l = l.cwiseMin(r).maxCoeff();
	t_r = l.cwiseMax(r).minCoeff();
	tMin = t_l;
	tMax = t_r;
	if (t_l > t_r + EPS)
		return false;
	if (t_r < -EPS)
		return false;
	return true;
}

std::pair<bool, float> Bounds3::IntersectWithTime(const Ray &ray) const {
	// invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because
	// Multiply is faster that Division
	float t_l = float_min, t_r = float_max;
	Vec3f l = (pMin - ray.origin).cwiseProduct(ray.inv_direction);
	Vec3f r = (pMax - ray.origin).cwiseProduct(ray.inv_direction);
	t_l = l.cwiseMin(r).maxCoeff();
	t_r = l.cwiseMax(r).minCoeff();
	if (t_l > t_r + EPS)
		return {false, float_max};
	if (t_r < -EPS)
		return {false, float_max};
	return {true, t_l};
}

Bounds3 Union(const Bounds3 &b1, const Bounds3 &b2) {
	return Bounds3{b1.pMin.cwiseMin(b2.pMin), b1.pMax.cwiseMax(b2.pMax)};
}

Bounds3 Union(const Bounds3 &b, const Vec3f &p) {
	return Bounds3{b.pMin.cwiseMin(p), b.pMax.cwiseMax(p)};
}