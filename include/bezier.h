//
// Created by ViXbob on 2022/12/18.
//

#ifndef CS171_ASSIGNMENT4_BEZIER_H
#define CS171_ASSIGNMENT4_BEZIER_H
#include <vector>

#include "aabb.h"
#include "core.h"
#include "geometry.h"
#include "kdtree_accel.h"

enum SamplingMode {
	Uniform,
	Adaptive,
};

class CurveSegment {
public:
	CurveSegment() = default;
	~CurveSegment() = default;

	CurveSegment(const std::vector<Vec2f> &_controlPoints,
	             const std::vector<float> &_weight,
	             const std::vector<float> &_knot, int _n, int _k, bool _rational, bool _reversed);
	[[nodiscard]] Bounds2 getBound();

	std::vector<Vec2f> controlPoints;
	std::vector<float> weights;
	std::vector<float> knots;
	int n{};
	int order{};
	bool rational{};
	bool reversed{};
};

class TrimList {
public:
	TrimList() = default;
	~TrimList() = default;

	explicit TrimList(std::shared_ptr<LoopedTrimCurve> &trimCurve);
	[[nodiscard]] Bounds2 getBound() const;

	void remove(std::shared_ptr<LoopedTrimCurve> &newTrim);
	void add(std::shared_ptr<LoopedTrimCurve> &newTrim);

	std::vector<std::shared_ptr<LoopedTrimCurve>> trims;
	Bounds2 bound;
	bool isClockwise{};
};


class LoopedTrimCurve {
public:
	LoopedTrimCurve() = default;
	~LoopedTrimCurve() = default;

	explicit LoopedTrimCurve(const std::vector<CurveSegment> &curve_segments);

	void findExtrema();
	void checkOrientation();
	[[nodiscard]] Bounds2 getBound() const;

	std::vector<CurveSegment> curveSegments;
	std::vector<Vec2f> controlPoints;
	std::vector<float> uExtrema, vExtrema;
	std::shared_ptr<TrimList> trimList;
	bool isClockwise{};
	Bounds2 bound;
};

class NURBS : public std::enable_shared_from_this<NURBS>, public Object {
public:
	NURBS() = default;
	~NURBS() override = default;

	NURBS(int m, int n, int _k, int _l);
	NURBS(int m, int n, int k, int l,
	      const std::vector<std::vector<Vec3f>> &_controlPoints,
	      const std::vector<std::vector<float>> &_weight,
	      const std::vector<float> &_knotM, const std::vector<float> &_knotN,
	      std::shared_ptr<BSDF> &bsdf);
	NURBS(int m, int n, int k, int l,
	      const std::vector<std::vector<Vec3f>> &_controlPoints,
	      const std::vector<float> &_knotM, const std::vector<float> &_knotN,
	      std::shared_ptr<BSDF> &bsdf);

	std::vector<std::vector<Vec3f>> controlPoints;
	std::vector<std::vector<float>> weight;
	std::vector<float> knotM, knotN;
	int k{}, l{}, u_m{}, v_m{}, u_n{}, v_n{}, u_p{}, v_p{};
	Bounds3 bound;
	std::shared_ptr<BSDF> bsdf;
	BVHAccelPtr bvh;
	KDTreeAccelPtr kdtree;
	std::vector<std::shared_ptr<IntervalObject>> interval_objects;
	std::vector<std::shared_ptr<LoopedTrimCurve>> looped_trim_curves;
	std::shared_ptr<TrimList> trimList{nullptr};

	void setControlPoint(int i, int j, Vec3f point);
	void setControlPoint(const std::vector<std::vector<Vec3f>> &_controlPoints);
	void setWeight(int i, int j, float w);
	void setWeight(const std::vector<std::vector<float>> &w);
	void setKnotM(const std::vector<float> &knot);
	void setKnotN(const std::vector<float> &knot);
	void setMaterial(std::shared_ptr<BSDF> &new_bsdf);
	void refineAndInitIntervalObject();
	void buildBVH();
	void buildKDTree();
	void preprocessTrimCurves();
	void init();
	void
	setTrimCurve(const std::vector<std::shared_ptr<LoopedTrimCurve>> &trimCurves);

	static std::pair<float, float> evaluateN(std::vector<float> &knot, float t,
	                                         int i, int k);
	Vertex evaluateWithNormal(float u, float v);
	Vertex evaluateWithNormal(double u, double v);
	Vertex evaluateWithNormalOld(float u, float v);
	void refine();
	std::vector<Vec4f> getHomoControlPoints(int i, int j) const;
	std::shared_ptr<TriangleMesh> generateMesh(SamplingMode mode = Uniform,
	                                           int sampleMSize = 100,
	                                           int sampleNSize = 100);
	Bounds3 getBounds() const override;
	bool intersect(const Ray &ray, Interaction &interaction) const override;
	static void intersectWithTrimCurve(const Vec2f &p, std::shared_ptr<TrimList> &tl, bool &keep);
};

class IntervalObject : public Object {
public:
	IntervalObject() = default;
	IntervalObject(std::shared_ptr<NURBS> &_surface, int _i, int _j,
	               std::shared_ptr<BSDF> bsdf = nullptr);
	std::shared_ptr<BSDF> bsdf;
	std::shared_ptr<NURBS> surface;
	Bounds3 bound;
	bool isTrimmed{};
	// surface patch [knotU[i], knotU[i + 1]) x [knotV[j], knotV[j + 1])
	int i{}, j{};
	void updateBounds();
	[[nodiscard]] Bounds3 getBounds() const override;
	bool intersect(const Ray &ray, Interaction &interaction) const override;
};

std::vector<std::vector<Vec3f>>
readControlPoints(const std::string &path, int m, int n,
                  const Vec3f &t = Vec3f(0.f, 0.f, 0.f));
std::vector<std::vector<float>> readWeights(const std::string &path, int m,
                                            int n);
void insertKnot(std::vector<float> &knot, std::vector<Vec4f> &controlPoints,
                int &p, int &m, int &n, float u);
Bounds3 getCurveBounds(const std::vector<float> &knot,
                       const std::vector<Vec4f> &controlPoints, const int &p,
                       const int &m, const int &n, int k);
void refineToBezier(std::vector<float> &knot, std::vector<Vec4f> &controlPoints,
                    int &p, int &m, int &n, int k, int base);

#endif// CS171_ASSIGNMENT4_BEZIER_H
