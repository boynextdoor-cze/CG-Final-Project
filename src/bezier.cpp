//
// Created by ViXbob on 2022/12/18.
//
#include "bezier.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>

#include "utils.h"

const float epsN = 1e-9;

// Helper functions
void testSurfacePatchBound(std::shared_ptr<NURBS> &surface, int i, int j);
Vec4f toHomogeneous(Vec3f p, float w);
Vec3f fromHomogeneous(Vec4f p);
Vec4d toHomogeneous(Vec3d p, double w);
Vec3d fromHomogeneous(Vec4d p);
Vec3f to3f(Vec3d p);
Vec3d to3d(Vec3f p);
float multiplyForBSplineBasisFunction(float a, float b, float c);

Vec4f toHomogeneous(Vec3f p, float w) {
  p *= w;
  return {p.x(), p.y(), p.z(), w};
}

Vec3f fromHomogeneous(Vec4f p) {
  p /= p.w();
  return p.head<3>();
}

Vec4d toHomogeneous(Vec3d p, double w) {
  p *= w;
  return {p.x(), p.y(), p.z(), w};
}

Vec3d fromHomogeneous(Vec4d p) {
  p /= p.w();
  return p.head<3>();
}

Vec3f to3f(Vec3d p) { return {(float)p.x(), (float)p.y(), (float)p.z()}; }

Vec3d to3d(Vec3f p) { return {(double)p.x(), (double)p.y(), (double)p.z()}; }

NURBS::NURBS(int m, int n, int _k, int _l) {
  controlPoints.resize(m + 1);
  for (auto &tmpControlPoints : controlPoints) {
    tmpControlPoints.resize(n + 1);
  }
  // initialize weight with 1.0f
  weight =
      std::vector<std::vector<float>>(m + 1, std::vector<float>(n + 1, 1.0f));
  k = _k;
  l = _l;
  knotM.resize(m + k + 1);
  knotN.resize(n + l + 1);
  u_m = m + k;
  u_n = m;
  u_p = _k - 1;

  v_m = n + l;
  v_n = n;
  v_p = _l - 1;
}

NURBS::NURBS(int _m, int _n, int _k, int _l,
             const std::vector<std::vector<Vec3f>> &_controlPoints,
             const std::vector<std::vector<float>> &_weight,
             const std::vector<float> &_knotM, const std::vector<float> &_knotN,
             std::shared_ptr<BSDF> &_bsdf)
    : u_m(_m + _k),
      u_n(_m),
      u_p(_k - 1),
      v_m(_n + _l),
      v_n(_n),
      v_p(_l - 1),
      k(_k),
      l(_l),
      controlPoints(_controlPoints),
      weight(_weight),
      knotM(_knotM),
      knotN(_knotN),
      bsdf(_bsdf) {
  if (u_m < 0 || u_n < 0 || u_p < 0 || v_m < 0 || v_n < 0 || v_p < 0 ||
      controlPoints.size() != u_n + 1 ||
      controlPoints.back().size() != v_n + 1 || weight.size() != u_n + 1 ||
      weight.back().size() != v_n + 1 || knotM.size() != u_m + 1 ||
      knotN.size() != v_m + 1) {
    std::cout << "ERROR::NURBS: Initialize error!" << std::endl;
    exit(-1);
  }
}

NURBS::NURBS(int _m, int _n, int _k, int _l,
             const std::vector<std::vector<Vec3f>> &_controlPoints,
             const std::vector<float> &_knotM, const std::vector<float> &_knotN,
             std::shared_ptr<BSDF> &_bsdf)
    : u_m(_m + _k),
      u_n(_m),
      u_p(_k - 1),
      v_m(_n + _l),
      v_n(_n),
      v_p(_l - 1),
      k(_k),
      l(_l),
      controlPoints(_controlPoints),
      knotM(_knotM),
      knotN(_knotN),
      bsdf(_bsdf) {
  weight = std::vector<std::vector<float>>(u_n + 1,
                                           std::vector<float>(v_n + 1, 1.0f));

  if (u_m < 0 || u_n < 0 || u_p < 0 || v_m < 0 || v_n < 0 || v_p < 0 ||
      controlPoints.size() != u_n + 1 ||
      controlPoints.back().size() != v_n + 1 || weight.size() != u_n + 1 ||
      weight.back().size() != v_n + 1 || knotM.size() != u_m + 1 ||
      knotN.size() != v_m + 1) {
    std::cout << "ERROR::NURBS: Initialize error!" << std::endl;
    exit(-1);
  }
}

void NURBS::setControlPoint(int i, int j, Vec3f point) {
  controlPoints[i][j] = point;
}

void NURBS::setControlPoint(
    const std::vector<std::vector<Vec3f>> &_controlPoints) {
  assert(_controlPoints.size() == controlPoints.size());
  assert(_controlPoints.front().size() == controlPoints.front().size());
  controlPoints = _controlPoints;
}

void NURBS::setWeight(int i, int j, float w) { weight[i][j] = w; }

void NURBS::setWeight(const std::vector<std::vector<float>> &w) {
  assert(w.size() == weight.size());
  assert(w.front().size() == weight.front().size());
  weight = w;
}

void NURBS::setKnotM(const std::vector<float> &knot) {
  assert(knot.size() == knotM.size());
  knotM = knot;
}

void NURBS::setKnotN(const std::vector<float> &knot) {
  assert(knot.size() == knotN.size());
  knotN = knot;
}

void NURBS::setMaterial(std::shared_ptr<BSDF> &new_bsdf) { bsdf = new_bsdf; }

void NURBS::refineAndInitIntervalObject() {
  refine();
  std::shared_ptr<NURBS> surface = shared_from_this();
  for (int i = u_p; i < u_m - u_p; i++) {
    for (int j = v_p; j < v_m - v_p; j++) {
      if (std::fabs(knotM[i + 1] - knotM[i]) < EQ_EPS ||
          std::fabs(knotN[j + 1] - knotN[j]) < EQ_EPS)
        continue;
      interval_objects.push_back(
          std::make_shared<IntervalObject>(surface, i, j, bsdf));
      bound = Union(bound, interval_objects.back()->getBounds());
    }
  }
}

void NURBS::buildBVH() {
  std::vector<ObjectPtr> objects(interval_objects.size());
  for (int i = 0; i < objects.size(); i++) objects[i] = interval_objects[i];
  bvh = std::make_shared<BVHAccel>(objects, 5);
}

void NURBS::buildKDTree() {
  std::vector<ObjectPtr> objects(interval_objects.size());
  for (int i = 0; i < objects.size(); i++) objects[i] = interval_objects[i];
	kdtree = std::make_shared<KDTreeAccel>(objects);
}

void NURBS::init() {
  refineAndInitIntervalObject();
  // buildBVH();
	buildKDTree();
}

std::pair<float, float> NURBS::evaluateN(std::vector<float> &knot, float t,
                                         int i, int k) {
  assert(!knot.empty());
  auto m = knot.size() - 1;
  assert(i + k <= m);
  auto tmp = std::vector<float>(m);
  for (auto index = 0; index < m; index++) {
    if (t >= knot[index] && t < knot[index + 1])
      tmp[index] = 1;
    else
      tmp[index] = 0;
  }
  for (auto iteration = 2; iteration < k; iteration++) {
    for (auto index = 0; index + iteration <= m; index++) {
      tmp[index] = multiplyForBSplineBasisFunction(
                       tmp[index], t - knot[index],
                       knot[index + iteration - 1] - knot[index]) +
                   multiplyForBSplineBasisFunction(
                       tmp[index + 1], knot[index + iteration] - t,
                       knot[index + iteration] - knot[index + 1]);
    }
  }
  auto value = multiplyForBSplineBasisFunction(tmp[i], t - knot[i],
                                               knot[i + k - 1] - knot[i]) +
               multiplyForBSplineBasisFunction(tmp[i + 1], knot[i + k] - t,
                                               knot[i + k] - knot[i + 1]);
  auto derivative =
      multiplyForBSplineBasisFunction(tmp[i], static_cast<float>(k) - 1,
                                      knot[i + k - 1] - knot[i]) -
      multiplyForBSplineBasisFunction(tmp[i + 1], static_cast<float>(k) - 1,
                                      knot[i + k] - knot[i + 1]);
  return {value, derivative};
}

Vertex NURBS::evaluateWithNormalOld(float u, float v) {
  assert(!controlPoints.empty());
  auto m = controlPoints.size(), n = controlPoints.front().size();
  Vec3f point(0.0f, 0.0f, 0.0f), tangentU(0.0f, 0.0f, 0.0f),
      tangentV(0.0f, 0.0f, 0.0f);
  float pointFrac = 0.0f, tangentUFrac = 0.0f, tangentVFrac = 0.0f;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      auto uM = evaluateN(knotM, u, i, k);
      auto vN = evaluateN(knotN, v, j, l);
      // calculate point at (u, v)
      auto tmp = weight[i][j] * uM.first * vN.first;
      point += controlPoints[i][j] * tmp;
      pointFrac += tmp;
      // d(u, v) / du
      tmp = weight[i][j] * uM.second * vN.first;
      tangentU += controlPoints[i][j] * tmp;
      tangentUFrac += tmp;
      // d(u, v) / dv
      tmp = weight[i][j] * uM.first * vN.second;
      tangentV += controlPoints[i][j] * tmp;
      tangentVFrac += tmp;
    }
  }
  assert(std::fabs(pointFrac) > epsN);
  auto p = point / pointFrac;
  auto U =
      (tangentU * pointFrac - tangentUFrac * point) / (pointFrac * pointFrac);
  auto V =
      (tangentV * pointFrac - tangentVFrac * point) / (pointFrac * pointFrac);
  auto normal = U.cross(V).normalized();
  return {p, normal, U, V};
}

Vertex NURBS::evaluateWithNormal(float u, float v) {
  // return evaluateWithNormalOld(u, v);
  // u \in [knotM[c], knotM[c + 1])
  // v \in [knotN[d], knotN[d + 1])
  int p = u_p, q = v_p;
  int c = upper_bound(knotM.begin(), knotM.end(), u) - knotM.begin() - 1;
  int d = upper_bound(knotN.begin(), knotN.end(), v) - knotN.begin() - 1;
  int s = 0, t = 0;
  for (int i = c; i >= 0 && std::fabs(knotM[i] - u) < EQ_EPS; i--) {
    s++;
  }
  for (int i = d; i >= 0 && std::fabs(knotN[i] - v) < EQ_EPS; i--) {
    t++;
  }

  std::vector<Vec4f> homo_points_on_u(p - s + 1);
  std::vector<Vec4f> homo_points_on_v(q - t + 1);

  Vec4f u_control_point, u_control_point_next;
  Vec4f v_control_point, v_control_point_next;
  // Compute along u direction
  for (int i = c - p; i <= c - s; i++) {
    for (int j = d - q; j <= d - t; j++) {
      homo_points_on_v[j - (d - q)] =
          toHomogeneous(controlPoints[i][j], weight[i][j]);
    }
    for (int r = 1; r <= q - t; r++) {
      for (int _i = d - t; _i >= d - q + r; _i--) {
        float alpha = (v - knotN[_i]) / (knotN[_i + q - r + 1] - knotN[_i]);
        int index = _i - (d - q);
        homo_points_on_v[index] = alpha * homo_points_on_v[index] +
                                  (1 - alpha) * homo_points_on_v[index - 1];
      }
    }
    homo_points_on_u[i - (c - p)] = homo_points_on_v[d - t - (d - q)];
  }

  for (int r = 1; r < p - s; r++) {
    for (int i = c - s; i >= c - p + r; i--) {
      float alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
      int index = i - (c - p);
      homo_points_on_u[index] = alpha * homo_points_on_u[index] +
                                (1 - alpha) * homo_points_on_u[index - 1];
    }
  }

  {
    int r = p - s, i = c - s;
    float alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
    int index = i - (c - p);
    u_control_point_next = homo_points_on_u[index];
    u_control_point = alpha * homo_points_on_u[index] +
                      (1 - alpha) * homo_points_on_u[index - 1];
  }

  // Compute along v direction
  for (int j = d - q; j <= d - t; j++) {
    for (int i = c - p; i <= c - s; i++) {
      homo_points_on_u[i - (c - p)] =
          toHomogeneous(controlPoints[i][j], weight[i][j]);
    }
    for (int r = 1; r <= p - s; r++) {
      for (int _i = c - s; _i >= c - p + r; _i--) {
        float alpha = (u - knotM[_i]) / (knotM[_i + p - r + 1] - knotM[_i]);
        int index = _i - (c - p);
        homo_points_on_u[index] = alpha * homo_points_on_u[index] +
                                  (1 - alpha) * homo_points_on_u[index - 1];
      }
    }
    homo_points_on_v[j - (d - q)] = homo_points_on_u[c - s - (c - p)];
  }

  for (int r = 1; r < q - t; r++) {
    for (int i = d - t; i >= d - q + r; i--) {
      float alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
      int index = i - (d - q);
      homo_points_on_v[index] = alpha * homo_points_on_v[index] +
                                (1 - alpha) * homo_points_on_v[index - 1];
    }
  }

  {
    int r = q - t, i = d - t;
    float alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
    int index = i - (d - q);
    v_control_point_next = homo_points_on_v[index];
    v_control_point = alpha * homo_points_on_v[index] +
                      (1 - alpha) * homo_points_on_v[index - 1];
  }

  Vec3f point = fromHomogeneous(u_control_point);
  Vec3f v_point = fromHomogeneous(v_control_point);
  Vec3f u_point_next = fromHomogeneous(u_control_point_next);
  Vec3f v_point_next = fromHomogeneous(v_control_point_next);
  Vec3f derivative_u = p * u_control_point_next.w() /
                       ((knotM[c + 1] - u) * u_control_point.w()) *
                       (u_point_next - point);
  Vec3f derivative_v = q * v_control_point_next.w() /
                       ((knotN[d + 1] - v) * u_control_point.w()) *
                       (v_point_next - point);
  Vec3f normal = derivative_u.cross(derivative_v).normalized();

// #define MY_DEBUG
#ifdef MY_DEBUG
  if ((point - v_point).norm() > EPS) {
    std::cout << "ERROR: evaluate self error!!" << std::endl;
    exit(-1);
  }
  Vertex evaluate_old = evaluateWithNormalOld(u, v);
  // float EPS = 1e-5;
  if ((evaluate_old.position - point).norm() > EPS) {
    exit(-1);
  }
  // std::cout << "Position" << std::endl;
  // std::cout << evaluate_old.position << std::endl;
  // std::cout << point << std::endl;
  // std::cout << "Derivative on u" << std::endl;
  // std::cout << evaluate_old.derivative_u << std::endl;
  // std::cout << derivative_u << std::endl;
  // std::cout << "Derivative on v" << std::endl;
  // std::cout << evaluate_old.derivative_v << std::endl;
  // std::cout << derivative_v << std::endl;
  // std::cout << "Normal" << std::endl;
  // std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v) <<
  // std::endl; std::cout << derivative_u.cross(derivative_v) << std::endl;
  // if((evaluate_old.derivative_u - derivative_u).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate derivative_u error!!" << std::endl;
  // 	exit(-1);
  // }
  // if((evaluate_old.derivative_v - derivative_v).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate derivate_v error!!" << std::endl;
  // 	std::cout << "Derivative on v" << std::endl;
  // 	std::cout << evaluate_old.derivative_v << std::endl;
  // 	std::cout << derivative_v << std::endl;
  // 	exit(-1);
  // }
  // if((evaluate_old.normal - normal).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate normal error!!" << std::endl;
  // 	std::cout << u << " " << v << std::endl;
  // 	// exit(-1);
  // }
  if ((evaluate_old.normal - normal).norm() > EPS ||
      (evaluate_old.derivative_v - derivative_v).norm() > EPS ||
      (evaluate_old.derivative_u - derivative_u).norm() > EPS) {
    std::cout << "ERROR: evaluate position error!!" << std::endl;
    std::cout << "Position" << std::endl;
    std::cout << evaluate_old.position << std::endl;
    std::cout << point << std::endl;
    std::cout << "Derivative on u" << std::endl;
    std::cout << evaluate_old.derivative_u << std::endl;
    std::cout << derivative_u << std::endl;
    std::cout << "Derivative on v" << std::endl;
    std::cout << evaluate_old.derivative_v << std::endl;
    std::cout << derivative_v << std::endl;
    std::cout << "Normal" << std::endl;
    std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v)
              << std::endl;
    std::cout << derivative_u.cross(derivative_v) << std::endl;
    std::cout << u << " " << v << std::endl;
    std::cout << s << " " << t << std::endl;
  }
#endif
  // #undef MY_DEBUG
  return {point, normal, derivative_u, derivative_v};
}

Vertex NURBS::evaluateWithNormal(double u, double v) {
  // return evaluateWithNormalOld(u, v);
  // u \in [knotM[c], knotM[c + 1])
  // v \in [knotN[d], knotN[d + 1])
  int p = u_p, q = v_p;
  int c = upper_bound(knotM.begin(), knotM.end(), u) - knotM.begin() - 1;
  int d = upper_bound(knotN.begin(), knotN.end(), v) - knotN.begin() - 1;
  int s = 0, t = 0;
  for (int i = c; i >= 0 && std::fabs(knotM[i] - u) < EQ_EPS; i--) {
    s++;
  }
  for (int i = d; i >= 0 && std::fabs(knotN[i] - v) < EQ_EPS; i--) {
    t++;
  }

  std::vector<Vec4d> homo_points_on_u(p - s + 1);
  std::vector<Vec4d> homo_points_on_v(q - t + 1);

  Vec4d u_control_point, u_control_point_next;
  Vec4d v_control_point, v_control_point_next;
  // Compute along u direction
  for (int i = c - p; i <= c - s; i++) {
    for (int j = d - q; j <= d - t; j++) {
      homo_points_on_v[j - (d - q)] =
          toHomogeneous(to3d(controlPoints[i][j]), (double)weight[i][j]);
    }
    for (int r = 1; r <= q - t; r++) {
      for (int _i = d - t; _i >= d - q + r; _i--) {
        double alpha = (v - knotN[_i]) / (knotN[_i + q - r + 1] - knotN[_i]);
        int index = _i - (d - q);
        homo_points_on_v[index] = alpha * homo_points_on_v[index] +
                                  (1 - alpha) * homo_points_on_v[index - 1];
      }
    }
    homo_points_on_u[i - (c - p)] = homo_points_on_v[d - t - (d - q)];
  }

  for (int r = 1; r < p - s; r++) {
    for (int i = c - s; i >= c - p + r; i--) {
      double alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
      int index = i - (c - p);
      homo_points_on_u[index] = alpha * homo_points_on_u[index] +
                                (1 - alpha) * homo_points_on_u[index - 1];
    }
  }

  {
    int r = p - s, i = c - s;
    double alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
    int index = i - (c - p);
    u_control_point_next = homo_points_on_u[index];
    u_control_point = alpha * homo_points_on_u[index] +
                      (1 - alpha) * homo_points_on_u[index - 1];
  }

  // Compute along v direction
  for (int j = d - q; j <= d - t; j++) {
    for (int i = c - p; i <= c - s; i++) {
      homo_points_on_u[i - (c - p)] =
          toHomogeneous(to3d(controlPoints[i][j]), (double)weight[i][j]);
    }
    for (int r = 1; r <= p - s; r++) {
      for (int _i = c - s; _i >= c - p + r; _i--) {
        double alpha = (u - knotM[_i]) / (knotM[_i + p - r + 1] - knotM[_i]);
        int index = _i - (c - p);
        homo_points_on_u[index] = alpha * homo_points_on_u[index] +
                                  (1 - alpha) * homo_points_on_u[index - 1];
      }
    }
    homo_points_on_v[j - (d - q)] = homo_points_on_u[c - s - (c - p)];
  }

  for (int r = 1; r < q - t; r++) {
    for (int i = d - t; i >= d - q + r; i--) {
      double alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
      int index = i - (d - q);
      homo_points_on_v[index] = alpha * homo_points_on_v[index] +
                                (1 - alpha) * homo_points_on_v[index - 1];
    }
  }

  {
    int r = q - t, i = d - t;
    double alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
    int index = i - (d - q);
    v_control_point_next = homo_points_on_v[index];
    v_control_point = alpha * homo_points_on_v[index] +
                      (1 - alpha) * homo_points_on_v[index - 1];
  }

  Vec3d point = fromHomogeneous(u_control_point);
  Vec3d v_point = fromHomogeneous(v_control_point);
  Vec3d u_point_next = fromHomogeneous(u_control_point_next);
  Vec3d v_point_next = fromHomogeneous(v_control_point_next);
  Vec3d derivative_u = p * u_control_point_next.w() /
                       ((knotM[c + 1] - u) * u_control_point.w()) *
                       (u_point_next - point);
  Vec3d derivative_v = q * v_control_point_next.w() /
                       ((knotN[d + 1] - v) * u_control_point.w()) *
                       (v_point_next - point);
  Vec3d normal = derivative_u.cross(derivative_v).normalized();

// #define MY_DEBUG
#ifdef MY_DEBUG
  if ((point - v_point).norm() > EPS) {
    std::cout << "ERROR: evaluate self error!!" << std::endl;
    exit(-1);
  }
  Vertex evaluate_old = evaluateWithNormalOld(u, v);
  // double EPS = 1e-5;
  if ((evaluate_old.position - to3f(point)).norm() > EPS) {
    exit(-1);
  }
  // std::cout << "Position" << std::endl;
  // std::cout << evaluate_old.position << std::endl;
  // std::cout << point << std::endl;
  // std::cout << "Derivative on u" << std::endl;
  // std::cout << evaluate_old.derivative_u << std::endl;
  // std::cout << derivative_u << std::endl;
  // std::cout << "Derivative on v" << std::endl;
  // std::cout << evaluate_old.derivative_v << std::endl;
  // std::cout << derivative_v << std::endl;
  // std::cout << "Normal" << std::endl;
  // std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v) <<
  // std::endl; std::cout << derivative_u.cross(derivative_v) << std::endl;
  // if((evaluate_old.derivative_u - derivative_u).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate derivative_u error!!" << std::endl;
  // 	exit(-1);
  // }
  // if((evaluate_old.derivative_v - derivative_v).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate derivate_v error!!" << std::endl;
  // 	std::cout << "Derivative on v" << std::endl;
  // 	std::cout << evaluate_old.derivative_v << std::endl;
  // 	std::cout << derivative_v << std::endl;
  // 	exit(-1);
  // }
  // if((evaluate_old.normal - normal).norm() > EPS) {
  // 	std::cout << "ERROR: evaluate normal error!!" << std::endl;
  // 	std::cout << u << " " << v << std::endl;
  // 	// exit(-1);
  // }
  float EPS = 1e-3;
  if ((evaluate_old.normal - to3f(normal)).norm() > EPS ||
      (evaluate_old.derivative_v - to3f(derivative_v)).norm() > EPS ||
      (evaluate_old.derivative_u - to3f(derivative_u)).norm() > EPS) {
    std::cout << "ERROR: evaluate position error!!" << std::endl;
    std::cout << "Position" << std::endl;
    std::cout << evaluate_old.position << std::endl;
    std::cout << point << std::endl;
    std::cout << "Derivative on u" << std::endl;
    std::cout << evaluate_old.derivative_u << std::endl;
    std::cout << derivative_u << std::endl;
    std::cout << "Derivative on v" << std::endl;
    std::cout << evaluate_old.derivative_v << std::endl;
    std::cout << derivative_v << std::endl;
    std::cout << "Normal" << std::endl;
    std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v)
              << std::endl;
    std::cout << derivative_u.cross(derivative_v) << std::endl;
    std::cout << u << " " << v << std::endl;
    std::cout << s << " " << t << std::endl;
    std::cout << (evaluate_old.normal - to3f(normal)).norm() << " "
              << (evaluate_old.derivative_v - to3f(derivative_v)).norm() << " "
              << (evaluate_old.derivative_u - to3f(derivative_u)).norm()
              << std::endl;
  }
#endif
  // #undef MY_DEBUG
  return {to3f(point), to3f(normal), to3f(derivative_u), to3f(derivative_v)};
}

void NURBS::refine() {
  static const int C = 35;
  std::vector<float> knotU_insert;
  std::vector<float> knotV_insert;
#ifdef MY_DEBUG
  std::cout << "start refinement" << std::endl;
  std::cout << "compute additional refined knots in each row" << std::endl;
#endif
  // calculate additional refined knots in each row
  {
    std::vector<Vec3f> V(u_n + 1), A(u_n + 1);
    for (int i = u_p; i < u_m - u_p; i++) {
      if (knotM[i] == knotM[i + 1]) continue;
      int max_n = 0;
      for (int ctrl_j = 0; ctrl_j <= v_n; ctrl_j++) {
        float max_A = 0;
        float sum_V = 0;
        for (int ctrl_i = i - u_p + 1; ctrl_i <= i; ctrl_i++)
          V[ctrl_i] = u_p / (knotM[ctrl_i + u_p] - knotM[ctrl_i]) *
                      (controlPoints[ctrl_i][ctrl_j] -
                       controlPoints[ctrl_i - 1][ctrl_j]);
        for (int ctrl_i = i - u_p + 2; ctrl_i <= i; ctrl_i++)
          A[ctrl_i] = (u_p - 1) / (knotM[ctrl_i + u_p - 1] - knotM[ctrl_i]) *
                      (V[ctrl_i] - V[ctrl_i - 1]);
        for (int j = i - u_p + 1; j <= i; j++) sum_V += V[j].norm();
        for (int j = i - u_p + 2; j <= i; j++)
          max_A = std::max(max_A, A[j].norm());
        sum_V = std::sqrt(sum_V / u_p);
        max_n = std::max(
            max_n,
            (int)std::ceil(C * max_A * std::pow(knotM[i + 1] - knotM[i], 1.5f) /
                           sum_V));
      }
      float segment = (knotM[i + 1] - knotM[i]) / (float)(max_n + 1);
      for (int j = 1; j <= max_n; j++)
        knotU_insert.push_back(knotM[i] + j * segment);
    }
  }
#ifdef MY_DEBUG
  std::cout << "compute additional refined knots in each column" << std::endl;
#endif
  // calculate additional refined knots in each column
  {
    std::vector<Vec3f> V(v_n + 1), A(v_n + 1);
    for (int i = v_p; i < v_m - v_p; i++) {
      if (knotN[i] == knotN[i + 1]) continue;
      int max_n = 0;
      for (int ctrl_i = 0; ctrl_i <= u_n; ctrl_i++) {
        float max_A = 0;
        float sum_V = 0;
        for (int ctrl_j = i - v_p + 1; ctrl_j <= i; ctrl_j++)
          V[ctrl_j] = v_p / (knotN[ctrl_j + v_p] - knotN[ctrl_j]) *
                      (controlPoints[ctrl_i][ctrl_j] -
                       controlPoints[ctrl_i][ctrl_j - 1]);
        for (int ctrl_j = i - v_p + 2; ctrl_j <= i; ctrl_j++)
          A[ctrl_j] = (v_p - 1) / (knotN[ctrl_j + v_p - 1] - knotN[ctrl_j]) *
                      (V[ctrl_j] - V[ctrl_j - 1]);
        for (int j = i - v_p + 1; j <= i; j++) sum_V += V[j].norm();
        for (int j = i - v_p + 2; j <= i; j++)
          max_A = std::max(max_A, A[j].norm());
        sum_V = std::sqrt(sum_V / v_p);
        max_n = std::max(
            max_n,
            (int)std::ceil(C * max_A * std::pow(knotN[i + 1] - knotN[i], 1.5f) /
                           sum_V));
      }
      float segment = (knotN[i + 1] - knotN[i]) / (float)(max_n + 1);
      for (int j = 1; j <= max_n; j++)
        knotV_insert.push_back(knotN[i] + j * segment);
    }
  }
#ifdef MY_DEBUG
  std::cout << "recompute control points while inserting knots" << std::endl;
#endif
  std::vector<std::vector<Vec4f>> homoControlPoints(
      u_n + 1, std::vector<Vec4f>(v_n + 1));

  for (int i = 0; i <= u_n; i++) {
    for (int j = 0; j <= v_n; j++)
      homoControlPoints[i][j] =
          toHomogeneous(controlPoints[i][j], weight[i][j]);
  }

  // Insert knots in knotV_insert and generate new control points
  for (int i = 0; i <= u_n; i++) {
    auto tmpKnot = knotN;
    auto tmpN = v_n;
    auto tmpM = v_m;
    for (float insert_knot : knotV_insert)
      insertKnot(tmpKnot, homoControlPoints[i], v_p, tmpM, tmpN, insert_knot);
    if (i == u_n) {
      knotN = tmpKnot;
      v_n = tmpN;
      v_m = tmpM;
    }
  }
#ifdef MY_DEBUG
  std::cout << "knots in knotV_insert have been inserted" << std::endl;
#endif
  // Insert knots in knotU_insert and generate new control points
  std::vector<std::vector<Vec4f>> newHomoControlPoints(
      v_n + 1, std::vector<Vec4f>(u_n + 1));
  for (int i = 0; i <= v_n; i++) {
    for (int j = 0; j <= u_n; j++) {
      newHomoControlPoints[i][j] = homoControlPoints[j][i];
    }
  }

  for (int i = 0; i <= v_n; i++) {
    auto tmpKnot = knotM;
    auto tmpN = u_n;
    auto tmpM = u_m;
    for (float insert_knot : knotU_insert)
      insertKnot(tmpKnot, newHomoControlPoints[i], u_p, tmpM, tmpN,
                 insert_knot);
    if (i == v_n) {
      knotM = tmpKnot;
      u_n = tmpN;
      u_m = tmpM;
    }
  }
#ifdef MY_DEBUG
  std::cout << "knots in knotU_insert have been inserted" << std::endl;
#endif
  // Ressign k and l
  k = u_p + 1;
  l = v_p + 1;

  // Reassign controlPoints
  controlPoints.resize(u_n + 1);
  weight.resize(u_n + 1);
  for (int i = 0; i <= u_n; i++) {
    controlPoints[i].resize(v_n + 1);
    weight[i].resize(v_n + 1);
    for (int j = 0; j <= v_n; j++) {
      controlPoints[i][j] = fromHomogeneous(newHomoControlPoints[j][i]);
      weight[i][j] = newHomoControlPoints[j][i].w();
    }
  }
// #define MY_DEBUG
#ifdef MY_DEBUG
  printf("u_m is %d, u_n is %d, u_p is %d\n", u_m, u_n, u_p);
  for (auto knot : knotM) std::cout << knot << " ";
  std::cout << std::endl;

  printf("v_m is %d, v_n is %d, v_p is %d\n", v_m, v_n, v_p);
  for (auto knot : knotN) std::cout << knot << " ";
  std::cout << std::endl;

  std::cout << "complete refinement" << std::endl;
#endif
// #undef MY_DEBUG

// #define MY_TEST_ON
#ifdef MY_TEST_ON
  {
    std::cout << "Test Start!" << std::endl;
    std::shared_ptr<NURBS> surface = shared_from_this();
    for (int i = u_p; i < u_m - u_p; i++) {
      for (int j = v_p; j < v_m - v_p; j++) {
        if (std::fabs(knotM[i + 1] - knotM[i]) < EQ_EPS ||
            std::fabs(knotN[j + 1] - knotN[j]) < EQ_EPS)
          continue;
        testSurfacePatchBound(surface, i, j);
      }
    }
  }
#endif
  // #undef MY_TEST_ON
}

std::vector<Vec4f> NURBS::getHomoControlPoints(int i, int j) const {
  if ((int)(i == -1) + (int)(j == -1) != 1) {
    std::cout
        << "ERROR::NURBS::getHomoControlPoints: Either i is -1 or j is -1!"
        << std::endl;
    exit(-1);
  }

  if (i == -1) {
    if (j < 0 || j > v_n) {
      std::cout << "ERROR::NURBS::getHomoControlPoints: j is out of bound!"
                << std::endl;
      exit(-1);
    }
    std::vector<Vec4f> result(u_n + 1);
    for (int _i = 0; _i <= u_n; _i++)
      result[_i] = toHomogeneous(controlPoints[_i][j], weight[_i][j]);
    return result;
  }

  if (j == -1) {
    if (i < 0 || i > u_n) {
      std::cout << "ERROR::NURBS::getHomoControlPoints: i is out of bound!"
                << std::endl;
      exit(-1);
    }
    std::vector<Vec4f> result(v_n + 1);
    for (int _j = 0; _j <= v_n; _j++)
      result[_j] = toHomogeneous(controlPoints[i][_j], weight[i][_j]);
    return result;
  }

  std::cout << "ERROR::NURBS::getHomoControlPoints: can not reach here!"
            << std::endl;
  exit(-1);
}

std::shared_ptr<TriangleMesh> NURBS::generateMesh(SamplingMode mode,
                                                  int sampleMSize,
                                                  int sampleNSize) {
  std::vector<Vec3f> vertices;
  std::vector<Vec3f> normals;
  std::vector<int> v_idx;
  if (mode == Uniform) {
    assert(sampleMSize > 1 && sampleNSize > 1);
    vertices.resize(sampleNSize * sampleMSize);
    normals.resize(sampleNSize * sampleMSize);
    float uL = knotM.front(), uR = knotM.back();
    float vL = knotN.front(), vR = knotN.back();
    for (int i = 0; i < sampleMSize; i++) {
      for (int j = 0; j < sampleNSize; j++) {
        auto u = uL + (uR - uL) * static_cast<float>(i) /
                          static_cast<float>(sampleMSize - 1);
        auto v = vL + (vR - vL) * static_cast<float>(j) /
                          static_cast<float>(sampleNSize - 1);
        // do some disturbance on corner sample
        if (i == 0) {
          u += 1e-5;
        }
        if (i == sampleMSize - 1) {
          u -= 1e-5;
        }
        if (j == 0) {
          v += 1e-5;
        }
        if (j == sampleNSize - 1) {
          v -= 1e-5;
        }
        auto vertex = evaluateWithNormal((double)u, (double)v);
        // vertices[i * sampleNSize + j] = 0.3f * vertex.position + Vec3f(-0.5,
        // 0.7f, 0);
        vertices[i * sampleNSize + j] = vertex.position;
        normals[i * sampleNSize + j] = -vertex.normal;
        // DEBUG_VEC(vertex.normal);
      }
    }

    for (int i = 0; i < sampleMSize - 1; i++) {
      for (int j = 0; j < sampleNSize - 1; j++) {
        v_idx.push_back(i * sampleNSize + j);
        v_idx.push_back(i * sampleNSize + j + 1);
        v_idx.push_back((i + 1) * sampleNSize + j);

        v_idx.push_back(i * sampleNSize + j + 1);
        v_idx.push_back((i + 1) * sampleNSize + j);
        v_idx.push_back((i + 1) * sampleNSize + j + 1);
      }
    }
  } else if (mode == Adaptive) {
    assert(false);
  } else {
    assert(false);
  }
  return std::make_shared<TriangleMesh>(vertices, normals, v_idx, v_idx);
}

Bounds3 NURBS::getBounds() const { return bound; }

bool NURBS::intersect(const Ray &ray, Interaction &interaction) const {
  if (bvh != nullptr) {
    // std::cout << "\rIntersect with NURBS" << std::endl;
    return bvh->getIntersection(ray, interaction);
  } 
	if (kdtree != nullptr) {
		return kdtree->getIntersection(ray, interaction);
	} else {
    if (!bound.IntersectP(ray)) return false;
    bool result = false;
    for (auto &interval_object : interval_objects) {
      result |= interval_object->intersect(ray, interaction);
    }
    return result;
  }
}

std::vector<std::vector<Vec3f>> readControlPoints(const std::string &path,
                                                  int m, int n, Vec3f translate) {
  std::vector<std::vector<Vec3f>> result(
      m + 1, std::vector<Vec3f>(n + 1, Vec3f(0.0f, 0.0f, 0.0f)));
  //  std::ifstream inputFile(getPath(path));
  std::ifstream inputFile(path);
  if (!inputFile.is_open()) {
    std::cout << "ERROR::Can not open the obj file!" << std::endl;
    return result;
  }
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      // inputFile >> result[i][j].x >> result[i][j].y >> result[i][j].z;
      inputFile >> result[i][j][0] >> result[i][j][1] >> result[i][j][2];
      result[i][j] = 0.1f * result[i][j] + translate;
    }
  }
  return std::move(result);
}

std::vector<std::vector<float>> readWeights(const std::string &path, int m,
                                            int n) {
  std::vector<std::vector<float>> result(m + 1,
                                         std::vector<float>(n + 1, 0.0f));
  // std::ifstream inputFile(getPath(path));
  std::ifstream inputFile(path);
  if (!inputFile.is_open()) {
    std::cout << "ERROR::Can not open the obj file!" << std::endl;
    return result;
  }
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      inputFile >> result[i][j];
    }
  }
  return std::move(result);
}

void insertKnot(std::vector<float> &knot, std::vector<Vec4f> &homoControlPoints,
                int &p, int &m, int &n, float u) {
  int k = upper_bound(knot.begin(), knot.end(), u) - knot.begin() - 1;
  int s = 0;
  for (int i = k; i >= 0 && std::fabs(knot[i] - u) < EQ_EPS; i--) s++;
  std::vector<Vec4f> q(p);
  for (int i = k - p + 1; i <= k - s; i++) {
    float alpha = (u - knot[i]) / (knot[i + p] - knot[i]);
    q[i - (k - p + 1)] =
        (1 - alpha) * homoControlPoints[i - 1] + alpha * homoControlPoints[i];
  }
  for (int i = k - p + 1; i < k - s; i++) {
    homoControlPoints[i] = q[i - (k - p + 1)];
  }
  homoControlPoints.insert(homoControlPoints.begin() + k - s,
                           q[k - (k - p + 1)]);
  knot.insert(knot.begin() + k + 1, u);
  m++;
  n++;
}

Bounds3 getCurveBounds(const std::vector<float> &knot,
                       const std::vector<Vec4f> &controlPoints, const int &p,
                       const int &m, const int &n, int k) {
  Bounds3 result;

  int c = k;
  int base = c - p;
  int s = 1;
  int h = p - s;
  float u = knot[k];

  std::vector<Vec4f> refinedControlPoints(controlPoints.begin() + c - p,
                                          controlPoints.begin() + c - s + 1);

  std::vector<Vec4f> q = refinedControlPoints;

  refinedControlPoints[p - 1] = controlPoints[k];
  for (int r = 1; r <= h; r++) {
    refinedControlPoints[p - 2 - (r - 1)] = q[c - s - base];
    for (int i = c - s; i >= c - p + r; i--) {
      float alpha = (u - knot[i]) / (knot[i + p - r + 1] - knot[i]);
      int index = i - base;
      q[index] = alpha * q[index] + (1 - alpha) * q[index - 1];
    }
  }
  result = Union(result, fromHomogeneous(q[c - s - base]));

  c = k + p;
  base = c - p;
  s = 1;
  h = p - s;
  u = knot[k + 1];
  q = refinedControlPoints;

  result = Union(result, fromHomogeneous(q[0]));
  for (int r = 1; r <= h; r++) {
    for (int i = c - s; i >= c - p + r; i--) {
      // knot[i]
      float knot_I = (i > k + p - 1) ? knot[k + i - (k + p - 1)] : knot[k];
      // knot[i + p - r + 1]
      float knot_J = (i + p - r + 1 > k + p - 1)
                         ? knot[k + (i + p - r + 1) - (k + p - 1)]
                         : knot[k];
      float alpha = (u - knot_I) / (knot_J - knot_I);
      int index = i - base;
      q[index] = alpha * q[index] + (1 - alpha) * q[index - 1];
    }
    result = Union(result, fromHomogeneous(q[c - p + r - base]));
  }

  return result;
}

IntervalObject::IntervalObject(std::shared_ptr<NURBS> &_surface, int _i, int _j,
                               std::shared_ptr<BSDF> _bsdf)
    : surface(_surface), i(_i), j(_j), bsdf(_bsdf) {
  updateBounds();
}

void refineToBezier(std::vector<float> &knot, std::vector<Vec4f> &controlPoints,
                    int &p, int &m, int &n, int k, int globalBase) {
  int c = k;
  int base = c - p;
  int s = 1;
  int h = p - s;
  float u = knot[k - globalBase];
  std::vector<Vec4f> refinedControlPoints(
      controlPoints.begin() + c - p - globalBase,
      controlPoints.begin() + c - s + 1 - globalBase);
  std::vector<Vec4f> q = refinedControlPoints;
  refinedControlPoints[p - 1] = controlPoints[k - globalBase];
  for (int r = 1; r <= h; r++) {
    refinedControlPoints[p - 2 - (r - 1)] = q[c - s - base];
    for (int i = c - s; i >= c - p + r; i--) {
      float alpha = (u - knot[i - globalBase]) /
                    (knot[i + p - r + 1 - globalBase] - knot[i - globalBase]);
      int index = i - base;
      q[index] = alpha * q[index] + (1 - alpha) * q[index - 1];
    }
    controlPoints[c - p + r - globalBase] = q[c - p + r - base];
  }

  c = k + p;
  base = c - p;
  s = 1;
  h = p - s;
  u = knot[k + 1 - globalBase];
  q = refinedControlPoints;

  // result = Union(result, fromHomogeneous(q[0]));
  refinedControlPoints[0] = q[0];
  for (int r = 1; r <= h; r++) {
    for (int i = c - s; i >= c - p + r; i--) {
      // knot[i]
      float knot_I = (i > k + p - 1) ? knot[k + i - (k + p - 1) - globalBase]
                                     : knot[k - globalBase];
      // knot[i + p - r + 1]
      float knot_J = (i + p - r + 1 > k + p - 1)
                         ? knot[k + (i + p - r + 1) - (k + p - 1) - globalBase]
                         : knot[k - globalBase];
      float alpha = (u - knot_I) / (knot_J - knot_I);
      int index = i - base;
      q[index] = alpha * q[index] + (1 - alpha) * q[index - 1];
    }
    refinedControlPoints[c - p + r - base] = q[c - p + r - base];
  }

  controlPoints.resize(2 * p);
  for (int i = k; i <= k + p - 1; i++) {
    controlPoints[i - globalBase] = refinedControlPoints[i - k];
  }
}

void IntervalObject::updateBounds() {
  bound = Bounds3();

  int uBase = i - surface->u_p, vBase = j - surface->v_p;

  std::vector<std::vector<Vec4f>> homoControlPoints(
      surface->u_p + 1, std::vector<Vec4f>(surface->v_p + 1));

  for (int index_i = i - surface->u_p; index_i <= i; index_i++) {
    for (int index_j = j - surface->v_p; index_j <= j; index_j++) {
      homoControlPoints[index_i - uBase][index_j - vBase] =
          toHomogeneous(surface->controlPoints[index_i][index_j],
                        surface->weight[index_i][index_j]);
    }
    std::vector<float> tmpKnotN(surface->knotN.begin() + j - surface->v_p,
                                surface->knotN.begin() + j + surface->v_p + 1);
    refineToBezier(tmpKnotN, homoControlPoints[index_i - uBase], surface->v_p,
                   surface->v_m, surface->v_n, j, vBase);
  }

  for (int index_j = j - 1; index_j <= j + surface->v_p - 1; index_j++) {
    std::vector<float> tmpKnotM(surface->knotM.begin() + i - surface->u_p,
                                surface->knotM.begin() + i + surface->u_p + 1);
    std::vector<Vec4f> finalHomoControlPoints(surface->u_p + 1);
    for (int index_i = i - surface->u_p; index_i <= i; index_i++)
      finalHomoControlPoints[index_i - uBase] =
          homoControlPoints[index_i - uBase][index_j - vBase];
    refineToBezier(tmpKnotM, finalHomoControlPoints, surface->u_p, surface->u_m,
                   surface->u_n, i, uBase);
    for (int index_i = i - 1; index_i <= i + surface->u_p - 1; index_i++) {
      bound = Union(bound,
                    fromHomogeneous(finalHomoControlPoints[index_i - uBase]));
    }
  }
}

Bounds3 IntervalObject::getBounds() const { return bound; }

struct NewtonIteration {
  NewtonIteration(double _u, double _v, const std::shared_ptr<NURBS> &surface,
                  const Ray &ray)
      : u(_u), v(_v) {
    p = surface->evaluateWithNormal(u, v);
    F_1 = ray.n_1.dot(p.position) + ray.d_1;
    F_2 = ray.n_2.dot(p.position) + ray.d_2;
    norm = std::sqrt(F_1 * F_1 + F_2 * F_2);
  }
  double u, v;
  double F_1, F_2;
  double norm;
  Vertex p;
};

bool IntervalObject::intersect(const Ray &ray, Interaction &interaction) const {
  static const int MAX_ITER = 7;
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> dis(0.0f, 1.0f);
  static const double EQ_EPS = 1e-10;
  static const double EPS = 1e-6;
  double u_0 = 0.5f * (surface->knotM[i] + surface->knotM[i + 1]);
  double v_0 = 0.5f * (surface->knotN[j] + surface->knotN[j + 1]);

  NewtonIteration curIteration(u_0, v_0, surface, ray);
  for (int iter = 0; iter < MAX_ITER; iter++) {
    if (curIteration.norm < EPS) {
      double t = (curIteration.p.position - ray.origin).dot(ray.direction);
      if (t < ray.t_min || t > ray.t_max) return false;
      if (t >= interaction.dist + EPS) return false;

      interaction.dist = t;
      interaction.pos = ray(t);
      interaction.normal = -curIteration.p.normal;
      interaction.material = bsdf;
      interaction.type = Interaction::Type::GEOMETRY;

      return true;
    }

    double J_11 = ray.n_1.dot(curIteration.p.derivative_u);
    double J_12 = ray.n_1.dot(curIteration.p.derivative_v);
    double J_21 = ray.n_2.dot(curIteration.p.derivative_u);
    double J_22 = ray.n_2.dot(curIteration.p.derivative_v);
    double det = J_11 * J_22 - J_12 * J_21;
    if (std::fabs(det) < EQ_EPS) {
      double next_u = curIteration.u + 0.1 * (u_0 - curIteration.u) * dis(gen);
      double next_v = curIteration.v + 0.1 * (v_0 - curIteration.v) * dis(gen);
      curIteration = NewtonIteration(next_u, next_u, surface, ray);
    } else {
      double invJ_11 = J_22, invJ_12 = -J_12, invJ_21 = -J_21, invJ_22 = J_11;
      double next_u =
          curIteration.u -
          (invJ_11 * curIteration.F_1 + invJ_12 * curIteration.F_2) / det;
      double next_v =
          curIteration.v -
          (invJ_21 * curIteration.F_1 + invJ_22 * curIteration.F_2) / det;

      if (next_u < surface->knotM[i] || next_u > surface->knotM[i + 1] ||
          next_v < surface->knotN[j] || next_v > surface->knotN[j + 1])
        return false;

      NewtonIteration nextIteration =
          NewtonIteration(next_u, next_v, surface, ray);
      if (nextIteration.norm - curIteration.norm > 0) return false;

      curIteration = nextIteration;
    }
  }
  return false;
}

// Test function
void testSurfacePatchBound(std::shared_ptr<NURBS> &surface, int i, int j) {
  std::cout << i << " " << j << std::endl;
  printf("[%f, %f) x [%f, %f)\n", surface->knotM[i], surface->knotM[i + 1],
         surface->knotN[j], surface->knotN[j + 1]);
  static const int TestBatch = 1000;
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<double> dis_u(surface->knotM[i],
                                               surface->knotM[i + 1]);
  std::uniform_real_distribution<double> dis_v(surface->knotN[j],
                                               surface->knotN[j + 1]);
  IntervalObject surface_patch(surface, i, j);
  auto bound = surface_patch.getBounds();
  for (int i = 0; i < TestBatch; i++) {
    // std::cout << "test " << i << std::endl;
    double u = dis_u(gen);
    double v = dis_v(gen);
    Vertex p = surface->evaluateWithNormal(u, v);
    if (!bound.Inside(p.position, bound)) {
      std::cout
          << "ERROR::Test_Failed: error bound, surface point is not in bound!"
          << std::endl;
      std::cout << u << " " << v << std::endl;
      std::cout << p.position << std::endl;
      std::cout << bound.pMin << std::endl;
      std::cout << bound.pMax << std::endl;
      exit(-1);
    }
  }
}

float multiplyForBSplineBasisFunction(float a, float b, float c) {
  if (std::fabs(c) < epsN) {
    assert(std::fabs(a) < epsN);
    return 0.0f;
  }
  return a * b / c;
}