#pragma once
#ifndef ACCEL_H_
#define ACCEL_H_

#include "bounds3.h"
#include "core.h"
#include "geometry.h"
#include "interaction.h"
#include "ray.h"
#include <vector>

// BVHAccel Forward Declarations
struct BVHObjectInfo;
struct LinearBVHNode {
	Bounds3 bounds;
	union {
		int object_offset;// leaf
		int right_offset; // interior
	};
	uint16_t objects_num;// 0 -> interior node
	uint8_t split_axis;  // interior node: xyz
	uint8_t padding;     // ensure 32 byte total size
};

class BVHAccel {

public:
	// BVHAccel Public Types
	enum class SplitMethod { NAIVE,
		                     SAH };

	// BVHAccel Public Methods
	explicit BVHAccel(std::vector<ObjectPtr> &p, int maxPrimsInNode = 1,
	                  SplitMethod splitMethod = SplitMethod::SAH);
	[[nodiscard]] Bounds3 WorldBound() const;
	~BVHAccel();

	bool getIntersection(const Ray &ray, Interaction &optimal,
	                     bool hitAny = false) const;
	void getIntersection(const BVHNodePtr &node, const Ray &ray,
	                     Interaction &optimal) const;
	BVHNodePtr root;

private:
	// BVHAccel Private Methods
	BVHNodePtr recursiveBuild(std::vector<BVHObjectInfo> &object_info, int begin,
	                          int end, int &total_nodes,
	                          std::vector<ObjectPtr> &ordered_obj);
	int flattenBVHTree(const BVHNodePtr &node, int &offset);
	void recursiveCheck(const BVHNodePtr &node);
	// BVHAccel Private Data
	const int maxPrimsInNode;
	const SplitMethod splitMethod;
	std::vector<ObjectPtr> objects;
	// std::vector<LinearBVHNode> nodes;
	LinearBVHNode *nodes;
};

struct BVHNode {
public:
	Bounds3 bounds;
	BVHNodePtr left;
	BVHNodePtr right;
	int split_axis, first_object_offset, objects_num;
	BVHNode();
	void initInterior(int _split_axis, BVHNodePtr _left, BVHNodePtr _right);
	void initLeaf(int first, int n, const Bounds3 &bound);
};

// You may need to add your code for BVH construction here.

#endif// ACCEL_H_