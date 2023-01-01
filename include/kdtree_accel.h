#pragma once
#ifndef KDTREE_ACCEL_H_
#define KDTREE_ACCEL_H_

#include "bounds3.h"
#include "core.h"
#include "geometry.h"
#include "interaction.h"
#include "ray.h"

struct KDAccelNode {
  // KdAccelNode Methods
  void InitLeaf(int *curObjectIndices, int np,
                std::vector<int> &objectInidces) {
    flags = 3;
    totalObjects |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
      objectIndex = 0;
    else if (np == 1)
      objectIndex = curObjectIndices[0];
    else {
      objectIndicesOffset = objectInidces.size();
      for (int i = 0; i < np; ++i) objectInidces.push_back(curObjectIndices[i]);
    }
  }
  void InitInterior(int axis, int ac, float s) {
    split = s;
    flags = axis;
    aboveChild |= (ac << 2);
  }
  float SplitPos() const { return split; }
  int TotalObjects() const { return totalObjects >> 2; }
  int SplitAxis() const { return flags & 3; }
  bool IsLeaf() const { return (flags & 3) == 3; }
  int AboveChild() const { return aboveChild >> 2; }
  // Bounds3 bound;
  union {
    float split;              // Interior
    int objectIndex;          // Leaf
    int objectIndicesOffset;  // Leaf
  };
  union {
    int flags;         // Both
    int totalObjects;  // Leaf
    int aboveChild;    // Interior
  };
};

enum class EdgeType { Start, End };
struct BoundEdge {
  // BoundEdge Public Methods
  BoundEdge() = default;
  BoundEdge(float t, int objectIndex, bool starting)
      : t(t), objectIndex(objectIndex) {
    type = starting ? EdgeType::Start : EdgeType::End;
  }
  float t;
  int objectIndex;
  EdgeType type;
};

class KDTreeAccel {
 public:
  KDTreeAccel() = default;
  ~KDTreeAccel() = default;
  KDTreeAccel(std::vector<ObjectPtr> &p, int intersectCost = 80,
              int traversalCost = 1, float emptyBonus = 0.5, int maxObjects = 2,
              int maxDepth = -1);
  Bounds3 getBounds() const { return bound; };
  bool getIntersection(const Ray &ray, Interaction &interaction) const;

 private:
  void recursiveBuild(const Bounds3 &curBound,
                      const std::vector<Bounds3> &objectBound,
                      int *curObjectIndices, int curTotalObjects, int depth,
                      std::vector<std::vector<BoundEdge>> &edges,
                      int *tmpObjects0, int *tmpObjects1, int badRefines);

  const int intersectCost, traversalCost, maxObjects;
  const float emptyBonus;
  std::vector<ObjectPtr> objects;
  std::vector<int> objectIndices;
  std::vector<KDAccelNode> nodes;
  Bounds3 bound;
};

struct KDToDo {
  const KDAccelNode *node;
  float tMin, tMax;
};
#endif  // KDTREE_ACCEL_H_