#pragma once
#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "aabb.h"
#include "accel.h"
#include "bsdf.h"
#include "core.h"
#include "interaction.h"
#include "ray.h"

#include <optional>
#include <vector>

struct Vertex {
	Vec3f position;
	Vec3f normal;
	Vec3f derivative_u, derivative_v;
	explicit Vertex(const Vec3f &_position = Vec3f(0, 0, 0),
	                const Vec3f &_normal = Vec3f(0, 0, 0),
	                const Vec3f &_derivative_u = Vec3f(0, 0, 0),
	                const Vec3f &_derivative_v = Vec3f(0, 0, 0)) {
		position = _position;
		normal = _normal;
		derivative_u = _derivative_u;
		derivative_v = _derivative_v;
	}
};

class Object {
public:
	Object() = default;
	virtual ~Object() = default;
	virtual bool intersect(const Ray &ray, Interaction &interaction) const = 0;
	[[nodiscard]] virtual Bounds3 getBounds() const = 0;
};

class Triangle : public Object {
public:
	Vec3f v0, v1, v2;// vertices A, B ,C , counter-clockwise order
	Vec3f n0, n1, n2;
	Bounds3 bound;
	std::shared_ptr<BSDF> bsdf;
	Triangle() = default;
	Triangle(const Vec3f &_v0, const Vec3f &_v1, const Vec3f &_v2, Vec3f _n0,
	         Vec3f _n1, Vec3f _n2, std::shared_ptr<BSDF> _bsdf = nullptr);
	Triangle(const std::vector<Vec3f> &vertices,
	         const std::vector<Vec3f> &normals, const Vec3i &v_idx,
	         const Vec3i &n_idx, std::shared_ptr<BSDF> _bsdf = nullptr);
	bool intersect(const Ray &ray, Interaction &interaction) const override;
	[[nodiscard]] Bounds3 getBounds() const override;
};

class TriangleMesh : public Object {
public:
	TriangleMesh() = default;
	TriangleMesh(std::vector<Vec3f> vertices, std::vector<Vec3f> normals,
	             std::vector<int> v_index, std::vector<int> n_index);
	bool intersect(const Ray &ray, Interaction &interaction) const override;
	[[nodiscard]] Bounds3 getBounds() const override;
	void setMaterial(std::shared_ptr<BSDF> &new_bsdf);
	void buildBVH();
	void pushAllTriangles(std::vector<ObjectPtr> &objects);

private:
	bool intersectOneTriangle(const Ray &ray, Interaction &interaction,
	                          const Vec3i &v_idx, const Vec3i &n_idx) const;
	std::shared_ptr<BSDF> bsdf;
	BVHAccelPtr bvh;
	std::vector<Vec3f> vertices;
	std::vector<Vec3f> normals;
	std::vector<int> v_indices;
	std::vector<int> n_indices;
	Bounds3 bound;
};

#endif// GEOMETRY_H_
