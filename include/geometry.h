#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "accel.h"
#include "bsdf.h"
#include "core.h"
#include "interaction.h"
#include "ray.h"

#include <optional>
#include <utility>
#include <vector>

struct Triangle {
	Triangle(unsigned int code,
	         std::vector<Vec3f> vs,
	         std::vector<Vec3f> ns,
	         std::shared_ptr<BSDF> material,
	         AABB aabb) : morton_code(code),
	                      vertices(std::move(vs)),
	                      normals(std::move(ns)),
	                      bsdf(std::move(material)),
	                      aabb(std::move(aabb)) {}
	unsigned int morton_code;
	std::shared_ptr<BSDF> bsdf;
	std::vector<Vec3f> vertices;
	std::vector<Vec3f> normals;
	AABB aabb;
	bool intersect(Ray &ray, Interaction &interaction);
};

class TriangleMesh {
public:
	TriangleMesh() = default;
	TriangleMesh(std::vector<Vec3f> vertices,
	             std::vector<Vec3f> normals,
	             std::vector<int> v_index,
	             std::vector<int> n_index);
	bool intersect(Ray &ray, Interaction &interaction) const;
	void setMaterial(std::shared_ptr<BSDF> &new_bsdf);
	AABB generateAABB();
	std::vector<std::shared_ptr<Triangle>> generateMortonCode(AABB &box);
	static unsigned int getMortonCode(const Vec3f &pos, const AABB &box);

private:
	bool intersectOneTriangle(Ray &ray, Interaction &interaction, const Vec3i &v_idx, const Vec3i &n_idx) const;
	std::shared_ptr<BSDF> bsdf;

	std::vector<Vec3f> vertices;
	std::vector<Vec3f> normals;
	std::vector<int> v_indices;
	std::vector<int> n_indices;
};

#endif// GEOMETRY_H_
