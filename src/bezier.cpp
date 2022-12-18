//
// Created by ViXbob on 2022/12/18.
//
#include "bezier.h"
#include "utils.h"
#include <fstream>
#include <iostream>

NURBS::NURBS(int m, int n, int _k, int _l) {
  controlPoints.resize(m + 1);
  for (auto &tmpControlPoints : controlPoints) {
    tmpControlPoints.resize(n + 1);
  }
  // initialize weight with 1.0f
  weight = std::vector<std::vector<float>>(m + 1, std::vector<float>(n + 1, 1.0f));
  k = _k;
  l = _l;
  knotM.resize(m + k + 1);
  knotN.resize(n + l + 1);
}

void NURBS::setControlPoint(int i, int j, Vec3f point) {
  controlPoints[i][j] = point;
}

void NURBS::setControlPoint(const std::vector<std::vector<Vec3f>> &_controlPoints) {
  assert(_controlPoints.size() == controlPoints.size());
  assert(_controlPoints.front().size() == controlPoints.front().size());
  controlPoints = _controlPoints;
}

void NURBS::setWeight(int i, int j, float w) {
  weight[i][j] = w;
}

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

const float epsN = 1e-9;

float multiplyForBSplineBasisFunction(float a, float b, float c) {
  if (std::fabs(c) < epsN) {
    assert(std::fabs(a) < epsN);
    return 0.0f;
  }
  return a * b / c;
}

std::pair<float, float> NURBS::evaluateN(std::vector<float> &knot, float t, int i, int k) {
  assert(!knot.empty());
  auto m = knot.size() - 1;
  assert(i + k <= m);
  auto tmp = std::vector<float>(m);
  for (auto index = 0; index < m; index++) {
    if (t >= knot[index] && t < knot[index + 1])
      tmp[index] = 1;
    else tmp[index] = 0;
  }
  for (auto iteration = 2; iteration < k; iteration++) {
    for (auto index = 0; index + iteration <= m; index++) {
      tmp[index] = multiplyForBSplineBasisFunction(tmp[index], t - knot[index], knot[index + iteration - 1] - knot[index])
                   + multiplyForBSplineBasisFunction(tmp[index + 1], knot[index + iteration] - t, knot[index + iteration] - knot[index + 1]);
    }
  }
  auto value = multiplyForBSplineBasisFunction(tmp[i], t - knot[i], knot[i + k - 1] - knot[i])
               + multiplyForBSplineBasisFunction(tmp[i + 1], knot[i + k] - t, knot[i + k] - knot[i + 1]);
  auto derivative = multiplyForBSplineBasisFunction(tmp[i], static_cast<float>(k) - 1, knot[i + k - 1] - knot[i])
                    - multiplyForBSplineBasisFunction(tmp[i + 1], static_cast<float>(k) - 1, knot[i + k] - knot[i + 1]);
  return {value, derivative};
}

Vertex NURBS::evaluateWithNormal(std::vector<std::vector<Vec3f>> &controlPoints,
                                                            std::vector<std::vector<float>> &weight,
                                                            std::vector<float> &knotM, std::vector<float> &knotN,
                                                            float u, float v, int k, int l) {
  assert(!controlPoints.empty());
  auto m = controlPoints.size(), n = controlPoints.front().size();
  Vec3f point(0.0f, 0.0f, 0.0f), tangentU(0.0f, 0.0f, 0.0f), tangentV(0.0f, 0.0f, 0.0f);
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
  auto U = (tangentU * pointFrac - tangentUFrac * point) / (pointFrac * pointFrac);
  auto V = (tangentV * pointFrac - tangentVFrac * point) / (pointFrac * pointFrac);
  auto normal = U.cross(V).normalized();
  return {p, normal};
}

std::shared_ptr<TriangleMesh> NURBS::generateMesh(SamplingMode mode, int sampleMSize, int sampleNSize) {
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
        auto u = uL + (uR - uL) * static_cast<float>(i) / static_cast<float>(sampleMSize - 1);
        auto v = vL + (vR - vL) * static_cast<float>(j) / static_cast<float>(sampleNSize - 1);
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
        auto vertex = NURBS::evaluateWithNormal(controlPoints, weight, knotM, knotN, u, v, k, l);
        vertices[i * sampleNSize + j] = 0.3f * vertex.position + Vec3f(-0.5, 0.7f, 0);
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

std::vector<std::vector<Vec3f>> readControlPoints(const std::string &path, int m, int n) {
  std::vector<std::vector<Vec3f>> result(m + 1, std::vector<Vec3f>(n + 1, Vec3f(0.0f, 0.0f, 0.0f)));
//  std::ifstream inputFile(getPath(path));
  std::ifstream inputFile(path);
  if(!inputFile.is_open()) {
    std::cout << "ERROR::Can not open the obj file!" << std::endl;
    return result;
  }
  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      // inputFile >> result[i][j].x >> result[i][j].y >> result[i][j].z;
      inputFile >> result[i][j][0] >> result[i][j][1] >> result[i][j][2];
    }
  }
  return std::move(result);
}

std::vector<std::vector<float>> readWeights(const std::string &path, int m, int n) {
  std::vector<std::vector<float>> result(m + 1, std::vector<float>(n + 1, 0.0f));
  // std::ifstream inputFile(getPath(path));
  std::ifstream inputFile(path);
  if(!inputFile.is_open()) {
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