#ifndef CORE_H_
#define CORE_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <memory>

#define UNIMPLEMENTED

template <typename T, int size>
using Vec = Eigen::Matrix<T, size, 1>;

using Vec3f = Eigen::Vector3f;
using Vec3i = Eigen::Vector3i;
using Vec2f = Eigen::Vector2f;
using Vec2i = Eigen::Vector2i;
using Vec4f = Eigen::Vector4f;
using Mat3f = Eigen::Matrix3f;
using Mat2d = Eigen::Matrix2d;
using Vec4d = Eigen::Vector4d;
using Vec3d = Eigen::Vector3d;
using Vec2d = Eigen::Vector2d;

constexpr float RAY_DEFAULT_MIN = 1e-5;
constexpr float RAY_DEFAULT_MAX = 1e7;
constexpr float PI = 3.141592653579f;
constexpr float INV_PI = 0.31830988618379067154;
constexpr float EPS = 1e-5;
constexpr float EQ_EPS = 1e-10;

class BSDF;
class Sampler;

class Bounds3;
class BVHNode;
class BVHAccel;
class Object;
class Vertex;
class NURBS;
class IntervalObject;
struct BVHObjectInfo;
struct LinearBVHNode;

using ObjectPtr = std::shared_ptr<Object>;
using BVHNodePtr = std::shared_ptr<BVHNode>;
using BVHAccelPtr = std::shared_ptr<BVHAccel>;

#endif //CORE_H_
