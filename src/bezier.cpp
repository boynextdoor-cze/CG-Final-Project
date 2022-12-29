//
// Created by ViXbob on 2022/12/18.
//
#include "bezier.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <algorithm>

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
	u_m = m + k;
	u_n = m;
	u_p = _k - 1;

	v_m = n + l;
	v_n = n;
	v_p = _l - 1;

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

Vertex NURBS::evaluateWithNormalOld(float u, float v) {
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
  return {p, normal, U, V};
}

Vertex NURBS::evaluateWithNormal(float u, float v) {
	return evaluateWithNormalOld(u, v);
	// u \in [knotM[c], knotM[c + 1])
	// v \in [knotN[d], knotN[d + 1])
	int p = u_p, q = v_p;
	int c = upper_bound(knotM.begin(), knotM.end(), u) - knotM.begin() - 1;
	int d = upper_bound(knotN.begin(), knotN.end(), v) - knotN.begin() - 1;
	int s = 0, t = 0;
	for(int i = c; i >= 0 && std::fabs(knotM[i] - u) < EQ_EPS; i--) {
		s++;
	}
	for(int i = d; i >= 0 && std::fabs(knotN[i] - v) < EQ_EPS; i--) {
		t++;
	}

	std::vector<Vec4f> homo_points_on_u(p - s + 1);
	std::vector<Vec4f> homo_points_on_v(q - t + 1);
	
	Vec4f u_control_point, u_control_point_next;
	Vec4f v_control_point, v_control_point_next;
	// Compute along u direction
	for(int i = c - p; i <= c - s; i++) {
		for(int j = d - q; j <= d - t; j++) {
			homo_points_on_v[j - (d - q)] = {controlPoints[i][j].x(), controlPoints[i][j].y(), controlPoints[i][j].z(), weight[i][j]};
		}
		for(int r = 1; r <= q - t; r++) {
			for(int _i = d - t; _i >= d - q + r; _i--) {
				float alpha = (v - knotN[_i]) / (knotN[_i + q - r + 1] - knotN[_i]);
				int index = _i - (d - q);
				homo_points_on_v[index] = alpha * homo_points_on_v[index] + (1 - alpha) * homo_points_on_v[index - 1];
			}
		}
		homo_points_on_u[i - (c - p)] = homo_points_on_v[d - t - (d - q)];
	}
	
	for(int r = 1; r < p - s; r++) {
		for(int i = c - s; i >= c - p + r; i--) {
			float alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
			int index = i - (c - p);
			homo_points_on_u[index] = alpha * homo_points_on_u[index] + (1 - alpha) * homo_points_on_u[index - 1];
		}
	}

	{
		int r = p - s, i = c - s;
		float alpha = (u - knotM[i]) / (knotM[i + p - r + 1] - knotM[i]);
		int index = i - (c - p);
		u_control_point_next = homo_points_on_u[index];
		u_control_point = alpha * homo_points_on_u[index] + (1 - alpha) * homo_points_on_u[index - 1];
	}

	// Compute along v direction
	for(int j = d - q; j <= d - t; j++) {
		for(int i = c - p; i <= c - s; i++) {
			homo_points_on_u[i - (c - p)] = {controlPoints[i][j].x(), controlPoints[i][j].y(), controlPoints[i][j].z(), weight[i][j]};
		}
		for(int r = 1; r <= p - s; r++) {
			for(int _i = c - s; _i >= c - p + r; _i--) {
				float alpha = (u - knotM[_i]) / (knotM[_i + p - r + 1] - knotM[_i]);
				int index = _i - (c - p);
				homo_points_on_u[index] = alpha * homo_points_on_u[index] + (1 - alpha) * homo_points_on_u[index - 1];
			}
		}
		homo_points_on_v[j - (d - q)] = homo_points_on_u[c - s - (c - p)];
	}

	for(int r = 1; r < q - t; r++) {
		for(int i = d - t; i >= d - q + r; i--) {
				float alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
				int index = i - (d - q);
				homo_points_on_v[index] = alpha * homo_points_on_v[index] + (1 - alpha) * homo_points_on_v[index - 1];
		}
	}

	{
		int r = q - t, i = d - t;
		float alpha = (v - knotN[i]) / (knotN[i + q - r + 1] - knotN[i]);
		int index = i - (d - q);
		v_control_point_next = homo_points_on_v[index];
		v_control_point = alpha * homo_points_on_v[index] + (1 - alpha) * homo_points_on_v[index - 1];
	}

	Vec3f point = u_control_point.head<3>() / u_control_point.w();
	Vec3f	v_point = v_control_point.head<3>() / v_control_point.w();
	Vec3f u_point_next = u_control_point_next.head<3>() / u_control_point_next.w();
	Vec3f v_point_next = v_control_point_next.head<3>() / v_control_point_next.w();
	Vec3f derivative_u = p * u_control_point_next.w() / ((knotM[c + 1] - u) * u_control_point.w()) * (u_point_next - point);
	Vec3f derivative_v = q * v_control_point_next.w() / ((knotN[d + 1] - v) * u_control_point.w()) * (v_point_next - point);
	Vec3f normal = derivative_u.cross(derivative_v).normalized();

// #define MY_DEBUG
#ifdef MY_DEBUG
	if((point - v_point).norm() > EPS) {
		std::cout << "ERROR: evaluate self error!!" << std::endl;
		exit(-1);
	}
	Vertex evaluate_old = evaluateWithNormalOld(controlPoints, weight, knotM, knotN, u, v, k, l);
	float EPS = 1e-3;
	if((evaluate_old.position - point).norm() > EPS) {
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
		std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v) << std::endl;
		std::cout << derivative_u.cross(derivative_v) << std::endl;
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
	// std::cout << evaluate_old.derivative_u.cross(evaluate_old.derivative_v) << std::endl;
	// std::cout << derivative_u.cross(derivative_v) << std::endl;
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
	if((evaluate_old.normal - normal).norm() > EPS || 
		 (evaluate_old.derivative_v - derivative_v).norm() > EPS || 
		 (evaluate_old.derivative_u - derivative_u).norm() > EPS) {
		std::cout << u << " " << v << std::endl;
		std::cout << s << " " << t << std::endl;
	}
#endif
// #undef MY_DEBUG
	return {point, normal, derivative_u, derivative_v};
}

void NURBS::refine() {
	static const int C = 3;
	std::vector<float> knotU_insert;
	std::vector<float> knotV_insert;
	
	// calculate additional refined knots in each row
	{
		std::vector<Vec3f> V(u_n + 1), A(u_n + 1);
		for(int i = u_p; i < u_m - u_p; i++) {
			if(knotM[i] == knotM[i + 1])
				continue;
			int max_n = 0;
			for(int ctrl_j = 0; ctrl_j <= v_n; ctrl_j++) {
				float max_A = 0;
				float sum_V = 0;
				for(int ctrl_i = i - u_p + 1; ctrl_i <= i; ctrl_i++) 
					V[ctrl_i] = u_p / (knotM[ctrl_i + u_p] - knotM[ctrl_i]) * (controlPoints[ctrl_i][ctrl_j] - controlPoints[ctrl_i - 1][ctrl_j]);
				for(int ctrl_i = i - u_p + 2; ctrl_i <= i; ctrl_i++) 
					A[ctrl_i] = (u_p - 1) / (knotM[ctrl_i + u_p - 1] - knotM[ctrl_i]) * (V[ctrl_i] - V[ctrl_i - 1]);
				for(int j = i - u_p + 1; j <= i; j++) 
					sum_V += V[j].norm();
				for(int j = i - u_p + 2; j <= i; j++) 
					max_A = std::max(max_A, A[j].norm());
				sum_V = std::sqrt(sum_V / u_p);
				max_n = std::max(max_n, (int)std::ceil(C * max_A * std::pow(knotM[i + 1] - knotM[i], 1.5f) / sum_V));
			}
			float segment = (knotM[i + 1] - knotM[i]) / (float)(max_n + 1);
			for(float new_knot = knotM[i] + segment; new_knot < knotM[i + 1]; new_knot += segment) {
				knotU_insert.push_back(new_knot);
			}
		}
	}

	// calculate additional refined knots in each column
	{
		std::vector<Vec3f> V(v_n + 1), A(v_n + 1);
		for(int i = v_p; i < v_m - v_p; i++) {
			if(knotN[i] == knotN[i + 1]) 
				continue;
			int max_n = 0;
			for(int ctrl_i = 0; ctrl_i <= u_n; ctrl_i++) {
				float max_A = 0;
				float sum_V = 0;
				for(int ctrl_j = i - v_p + 1; ctrl_j <= i; ctrl_j++)
					V[ctrl_j] = v_p / (knotN[ctrl_j + v_p] - knotN[ctrl_j]) * (controlPoints[ctrl_i][ctrl_j] - controlPoints[ctrl_i][ctrl_j - 1]);
				for(int ctrl_j = i - v_p + 2; ctrl_j <= i; ctrl_j++)
					A[ctrl_j] = (v_p - 1) / (knotN[ctrl_j + v_p - 1] - knotN[ctrl_j]) * (V[ctrl_j] - V[ctrl_j - 1]);
				for(int j = i - v_p + 1; j <= i; j++)
					sum_V += V[j].norm();
				for(int j = i - v_p + 2; j <= i; j++)
					max_A = std::max(max_A, A[j].norm());
				sum_V = std::sqrt(sum_V / v_p);
				max_n = std::max(max_n, (int)std::ceil(C * max_A * std::pow(knotN[i + 1] - knotN[i], 1.5f) / sum_V));
			}
			float segment = (knotN[i + 1] - knotN[i]) / (float)(max_n + 1);
			for(float new_knot = knotN[i] + segment; new_knot < knotN[i + 1]; new_knot += segment)
				knotV_insert.push_back(new_knot);
		}
	}

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
        auto vertex = evaluateWithNormal(u, v);
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