//
// Created by ViXbob on 2022/12/18.
//

#ifndef CS171_ASSIGNMENT4_BEZIER_H
#define CS171_ASSIGNMENT4_BEZIER_H
#include "core.h"
#include "geometry.h"
#include <vector>

enum SamplingMode {
    Uniform,
    Adaptive,
};

class NURBS {
public:
    std::vector<std::vector<Vec3f>> controlPoints;
    std::vector<std::vector<float>> weight;
    std::vector<float> knotM, knotN;
    int k, l, u_m, v_m, u_n, v_n, u_p, v_p;
    NURBS(int m, int n, int _k, int _l);
    void setControlPoint(int i, int j, Vec3f point);
    void setControlPoint(const std::vector<std::vector<Vec3f>> &_controlPoints);
    void setWeight(int i, int j, float w);
    void setWeight(const std::vector<std::vector<float>> &w);
    void setKnotM(const std::vector<float> &knot);
    void setKnotN(const std::vector<float> &knot);
    std::pair<float, float> evaluateN(std::vector<float> &knot, float t, int i, int k);
    Vertex evaluateWithNormal(float u, float v);
		Vertex evaluateWithNormal(double u, double v);
		Vertex evaluateWithNormalOld(float u, float v);
		void refine();
    std::shared_ptr<TriangleMesh> generateMesh(SamplingMode mode = Uniform, int sampleMSize = 100, int sampleNSize = 100);
};

std::vector<std::vector<Vec3f>> readControlPoints(const std::string &path, int m, int n);
std::vector<std::vector<float>> readWeights(const std::string &path, int m, int n);
void insertKnot(std::vector<float> &knot, std::vector<Vec4f> &controlPoints, int &p, int &m, int &n, float u);
#endif //CS171_ASSIGNMENT4_BEZIER_H
