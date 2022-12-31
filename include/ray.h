#ifndef RAY_H_
#define RAY_H_

#include "core.h"

struct Ray {
  /// origin point of ray
  Vec3f origin;
  /// normalized direction of the ray
  Vec3f direction;
  Vec3f inv_direction;
  /// min and max distance of the ray
  float t_min;
  float t_max;
	/// parameter for ray nurbs
	Vec3f n_1, n_2;
	float d_1, d_2;

  explicit Ray(const Vec3f &o, const Vec3f &dir, float t_min = RAY_DEFAULT_MIN, float t_max = RAY_DEFAULT_MAX)
  : origin(o)
  , direction(dir.normalized())
  , t_min(t_min)
  , t_max(t_max) {
		direction.normalize();
    inv_direction = {
      std::fabs(direction[0]) < (float)1e-10 ? (float)1e10 : 1.f / direction[0],
      std::fabs(direction[1]) < (float)1e-10 ? (float)1e10 : 1.f / direction[1],
      std::fabs(direction[2]) < (float)1e-10 ? (float)1e10 : 1.f / direction[2]
    };
		if(std::fabs(direction.x()) > std::fabs(direction.y()) && 
			 std::fabs(direction.x()) > std::fabs(direction.z())) {
			n_1 = {direction.y(), -direction.x(), 0.f};
			n_1.normalize();
		} else {
			n_1 = {0.f, direction.z(), -direction.y()};
			n_1.normalize();
		}
		n_2 = n_1.cross(direction).normalized();
		d_1 = -n_1.dot(origin);
		d_2 = -n_2.dot(origin);
  }

  [[nodiscard]] Vec3f operator()(float t) const {
    return origin + t * direction;
  }
};

#endif //RAY_H_
