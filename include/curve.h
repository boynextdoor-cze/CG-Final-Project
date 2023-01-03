//
// Created by boynextdoor on 1/3/23.
//

#ifndef CS171_FINAL_PROJECT_CURVE_H
#define CS171_FINAL_PROJECT_CURVE_H
#include "core.h"
#include <vector>

class TrimCurve {
public:
	TrimCurve() = default;
	~TrimCurve() = default;

	TrimCurve(const std::vector<Vec2f> &_controlPoints, const std::vector<float> &_weight, const std::vector<float> &_knot, int _n, int _k);

	std::vector<Vec2f> controlPoints;
	std::vector<float> weights;
	std::vector<float> knots;
	int n{};
	int order{};
};

#endif//CS171_FINAL_PROJECT_CURVE_H
