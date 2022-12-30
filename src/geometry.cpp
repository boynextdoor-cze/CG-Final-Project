#include "geometry.h"
#include "utils.h"

#include <cmath>
#include <iostream>
#include <utility>

bool isSigularTriangle(Vec3f v0, Vec3f v1, Vec3f v2) {
  if ((v1 - v0).cross(v2 - v0).norm() * 0.5f < 1e-15) {
    return true;
  } else
    return false;
}

Triangle::Triangle(Vec3f _v0, Vec3f _v1, Vec3f _v2, Vec3f _n0, Vec3f _n1,
                   Vec3f _n2, std::shared_ptr<BSDF> _bsdf)
    : v0(_v0), v1(_v1), v2(_v2), n0(_n0), n1(_n1), n2(_n2), bsdf(_bsdf),
      bound(_v0, _v1, _v2) {}

Triangle::Triangle(const std::vector<Vec3f> &vertices,
                   const std::vector<Vec3f> &normals, const Vec3i &v_idx,
                   const Vec3i &n_idx, std::shared_ptr<BSDF> _bsdf)
    : v0(vertices[v_idx[0]]), v1(vertices[v_idx[1]]), v2(vertices[v_idx[2]]),
      n0(normals[n_idx[0]]), n1(normals[n_idx[1]]), n2(normals[n_idx[2]]),
      bsdf(_bsdf),
      bound(vertices[v_idx[0]], vertices[v_idx[1]], vertices[v_idx[2]]) {}

bool Triangle::intersect(const Ray &ray, Interaction &interaction) const {
  Vec3f v0v1 = v1 - v0;
  Vec3f v0v2 = v2 - v0;
  Vec3f pvec = ray.direction.cross(v0v2);
  float det = v0v1.dot(pvec);

  float invDet = 1.0f / det;

  Vec3f tvec = ray.origin - v0;
  float u = tvec.dot(pvec) * invDet;
  if (u < 0 || u > 1)
    return false;
  Vec3f qvec = tvec.cross(v0v1);
  float v = ray.direction.dot(qvec) * invDet;
  if (v < 0 || u + v > 1)
    return false;
  float t = v0v2.dot(qvec) * invDet;
  if (t < ray.t_min || t > ray.t_max)
    return false;
  if (t >= interaction.dist + EPS)
    return false;

  auto normal = (u * n1 + v * n2 + (1 - u - v) * n0).normalized();

  if (std::isnan(normal.x()) ||
      std::isnan(normal.y()) ||
      std::isnan(normal.z())) {
    return false;
    DEBUG_VEC(ray.origin);
    DEBUG_VEC(ray.direction);
    DEBUG_VEC(v0);
    DEBUG_VEC(v1);
    DEBUG_VEC(v2);
  }

  interaction.dist = t;
  interaction.pos = ray(t);
  interaction.normal = normal;
  interaction.material = bsdf;
  interaction.type = Interaction::Type::GEOMETRY;
  return true;
}

Bounds3 Triangle::getBounds() const { return bound; }

TriangleMesh::TriangleMesh(std::vector<Vec3f> _vertices,
                           std::vector<Vec3f> _normals,
                           std::vector<int> _v_index, std::vector<int> _n_index)
    : vertices(std::move(_vertices)), normals(std::move(_normals)),
      v_indices(std::move(_v_index)), n_indices(std::move(_n_index)),
      bvh(nullptr) {
  for (int i = 0; i < v_indices.size() / 3; i++) {
    Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
    bound = Union(bound, Bounds3(vertices[v_idx[0]], vertices[v_idx[1]],
                                 vertices[v_idx[2]]));
  }
}

bool TriangleMesh::intersect(const Ray &ray, Interaction &interaction) const {
  if (bvh != nullptr) {
    // bvh->getIntersection(bvh->root, ray, interaction);
    return bvh->getIntersection(ray, interaction);
  } else {
    // If you did not implement BVH
    // directly loop through all triangles in the mesh and test intersection for
    // each triangle.
    for (int i = 0; i < v_indices.size() / 3; i++) {
      Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
      Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
      Interaction temp;
      if (intersectOneTriangle(ray, temp, v_idx, n_idx) &&
          (temp.dist < interaction.dist)) {
        interaction = temp;
      }
    }
  }
  return interaction.type != Interaction::Type::NONE;
}

Bounds3 TriangleMesh::getBounds() const { return bound; }

void TriangleMesh::setMaterial(std::shared_ptr<BSDF> &new_bsdf) {
  bsdf = new_bsdf;
}
void TriangleMesh::buildBVH() {
  assert(bsdf != nullptr);
  std::vector<ObjectPtr> objects(v_indices.size() / 3);
  int cnt = 0;
  for (int i = 0; i < v_indices.size() / 3; i++) {
    Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
    Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
    if (!isSigularTriangle(vertices[v_idx[0]], vertices[v_idx[1]],
                           vertices[v_idx[2]]))
      objects[cnt++] =
          std::make_shared<Triangle>(vertices, normals, v_idx, n_idx, bsdf);
  }
  objects.resize(cnt);
  bvh = std::make_shared<BVHAccel>(objects, 5);
}

void TriangleMesh::pushAllTriangles(std::vector<ObjectPtr> &objects) {
  for (int i = 0; i < v_indices.size() / 3; i++) {
    Vec3i v_idx(v_indices[3 * i], v_indices[3 * i + 1], v_indices[3 * i + 2]);
    Vec3i n_idx(n_indices[3 * i], n_indices[3 * i + 1], n_indices[3 * i + 2]);
    objects.push_back(
        std::make_shared<Triangle>(vertices, normals, v_idx, n_idx, bsdf));
  }
}

bool TriangleMesh::intersectOneTriangle(const Ray &ray,
                                        Interaction &interaction,
                                        const Vec3i &v_idx,
                                        const Vec3i &n_idx) const {
  Vec3f v0 = vertices[v_idx[0]];
  Vec3f v1 = vertices[v_idx[1]];
  Vec3f v2 = vertices[v_idx[2]];
  Vec3f v0v1 = v1 - v0;
  Vec3f v0v2 = v2 - v0;
  Vec3f pvec = ray.direction.cross(v0v2);
  float det = v0v1.dot(pvec);

  if (std::fabs(det) < 1e-15) {
    return false;
  }

  float invDet = 1.0f / det;

  Vec3f tvec = ray.origin - v0;
  float u = tvec.dot(pvec) * invDet;
  if (u < 0 || u > 1)
    return false;
  Vec3f qvec = tvec.cross(v0v1);
  float v = ray.direction.dot(qvec) * invDet;
  if (v < 0 || u + v > 1)
    return false;
  float t = v0v2.dot(qvec) * invDet;
  if (t < ray.t_min || t > ray.t_max)
    return false;

  interaction.dist = t;
  interaction.pos = ray(t);
  interaction.normal = (u * normals[n_idx[1]] + v * normals[n_idx[2]] +
                        (1 - u - v) * normals[n_idx[0]])
                           .normalized();
  interaction.material = bsdf;
  interaction.type = Interaction::Type::GEOMETRY;
  return true;
}