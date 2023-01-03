//
// Created by boynextdoor on 1/3/23.
//

#include "curve.h"

TrimCurve::TrimCurve(const std::vector<Vec2f> &_controlPoints, const std::vector<float> &_weight, const std::vector<float> &_knot, int _n, int _k) {
	controlPoints = _controlPoints;
	weights = _weight;
	knots = _knot;
	n = _n;
	order = _k;
}
