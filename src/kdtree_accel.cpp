#include "kdtree_accel.h"
#include <ctime>

KDTreeAccel::KDTreeAccel(std::vector<ObjectPtr> &_objects, int _intersectCost,
                         int _traversalCost, float _emptyBonus, int _maxObjects,
                         int maxDepth)
    : intersectCost(_intersectCost),
      traversalCost(_traversalCost),
      maxObjects(_maxObjects),
      emptyBonus(_emptyBonus),
      objects(_objects) {
  if (maxDepth <= 0)
    maxDepth = std::round(8 + 1.3f * std::log2(objects.size()));
  std::vector<Bounds3> objectBounds;
  objectBounds.reserve(objects.size());
  for (auto &object : objects) {
    auto curBound = object->getBounds();
    objectBounds.push_back(curBound);
    bound = Union(bound, curBound);
  }

  std::vector<std::vector<BoundEdge>> edges(
      3, std::vector<BoundEdge>(2 * objects.size()));
  std::vector<int> tmpObjects0(objects.size()),
      tmpObjects1((maxDepth + 1) * objects.size()), objectIndex(objects.size());
  for (int i = 0; i < objectIndex.size(); i++) objectIndex[i] = i;
  time_t start, stop;
  time(&start);
	printf("KDTree building starts\n");
  recursiveBuild(bound, objectBounds, objectIndex.data(), objects.size(),
                 maxDepth, edges, tmpObjects0.data(), tmpObjects1.data(), 0);
  time(&stop);
  double diff = difftime(stop, start);
  int hrs = (int)diff / 3600;
  int mins = ((int)diff / 60) - (hrs * 60);
  int secs = (int)diff - (hrs * 3600) - (mins * 60);
  printf(
      "\rKDTree Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
      hrs, mins, secs);
}

void KDTreeAccel::recursiveBuild(const Bounds3 &bound,
                                 const std::vector<Bounds3> &objectBound,
                                 int *curObjectIndices, int curTotalObjects,
                                 int depth,
                                 std::vector<std::vector<BoundEdge>> &edges,
                                 int *tmpObjects0, int *tmpObjects1,
                                 int badRefines) {
  if (curTotalObjects <= maxObjects || depth == 0) {
    nodes.emplace_back();
    nodes.back().InitLeaf(curObjectIndices, curTotalObjects, objectIndices);
    return;
  }
  int bestAxis = -1, bestOffset = -1;
  float bestCost = 1e9;
  float oldCost = intersectCost * (float)curTotalObjects;
  float totalSA = bound.SurfaceArea();
  float invTotalSA = 1.0f / totalSA;
  Vec3f d = bound.Diagonal();
  int axis = bound.maxExtent();
  int retries = 0;

retrySplit:
  for (int i = 0; i < curTotalObjects; i++) {
    int index = curObjectIndices[i];
    const Bounds3 &curBound = objectBound[index];
    edges[axis][2 * i] = BoundEdge(curBound.pMin[axis], index, true);
    edges[axis][2 * i + 1] = BoundEdge(curBound.pMax[axis], index, false);
  }

  std::sort(&edges[axis][0], &edges[axis][2 * curTotalObjects],
            [](const BoundEdge &edge_0, const BoundEdge &edge_1) -> bool {
              if (edge_0.t != edge_0.t)
                return edge_0.t < edge_1.t;
              else
                return (int)edge_0.type < (int)edge_1.type;
            });

  int totalBelow = 0, totalAbove = curTotalObjects;
  for (int i = 0; i < 2 * curTotalObjects; i++) {
    if (edges[axis][i].type == EdgeType::End) totalAbove--;
    float curEdgeT = edges[axis][i].t;

    if (curEdgeT > bound.pMin[axis] && curEdgeT < bound.pMax[axis]) {
      int otherAxis_0 = (axis + 1) % 3, otherAxis_1 = (axis + 2) % 3;
      float belowSA = d[otherAxis_0] * d[otherAxis_1] +
                      2.0f * (curEdgeT - bound.pMin[axis]) *
                          (d[otherAxis_0] + d[otherAxis_1]);
      float aboveSA = d[otherAxis_0] * d[otherAxis_1] +
                      2.0f * (bound.pMax[axis] - curEdgeT) *
                          (d[otherAxis_0] + d[otherAxis_1]);
      float pBelow = belowSA * invTotalSA;
      float pAbove = aboveSA * invTotalSA;
      float curEmptyBonus =
          (totalBelow == 0 || totalAbove == 0) ? emptyBonus : 0;
      float cost =
          traversalCost + intersectCost * (1 - emptyBonus) *
                              (pBelow * totalBelow + pAbove * totalAbove);

      if (cost < bestCost) {
        bestCost = cost;
        bestAxis = axis;
        bestOffset = i;
      }
      if (edges[axis][i].type == EdgeType::Start) totalBelow++;
    }
  }

  if (bestAxis == -1 && retries < 2) {
    ++retries;
    axis = (axis + 1) % 3;
    goto retrySplit;
  }

  if (bestCost > oldCost) ++badRefines;
  if ((bestCost > 4 * oldCost && curTotalObjects < 16) || bestAxis == -1 ||
      badRefines == 3) {
    nodes.emplace_back();
    nodes.back().InitLeaf(curObjectIndices, curTotalObjects, objectIndices);
    return;
  }

  int n0 = 0, n1 = 0;
  for (int i = 0; i < bestOffset; i++)
    if (edges[bestAxis][i].type == EdgeType::Start)
      tmpObjects0[n0++] = edges[bestAxis][i].objectIndex;
  for (int i = bestOffset + 1; i < 2 * curTotalObjects; i++)
    if (edges[bestAxis][i].type == EdgeType::End)
      tmpObjects1[n1++] = edges[bestAxis][i].objectIndex;

  float tSplit = edges[bestAxis][bestOffset].t;
  Bounds3 bound0 = bound, bound1 = bound;
  bound0.pMax[bestAxis] = bound1.pMin[bestAxis] = tSplit;
  nodes.emplace_back();
  auto &curNode = nodes.back();
  recursiveBuild(bound0, objectBound, tmpObjects0, n0, depth - 1, edges,
                 curObjectIndices, tmpObjects1 + n1, badRefines);
  curNode.InitInterior(bestAxis, nodes.size(), tSplit);
  recursiveBuild(bound1, objectBound, tmpObjects1, n1, depth - 1, edges,
                 curObjectIndices, tmpObjects1 + n1, badRefines);
}

bool KDTreeAccel::getIntersection(const Ray &ray,
                                  Interaction &interaction) const {
  float tMin, tMax;
  if (!bound.IntersectP(ray, tMin, tMax)) return false;

  KDToDo todo[64];
  int todoPos = 0;
  bool hit = false;
  const KDAccelNode *node = &nodes.data()[0];
  while (node != nullptr) {
    if (ray.t_max < tMin) break;

    if (!node->IsLeaf()) {
      int axis = node->SplitAxis();
      float tPlane =
          (node->SplitPos() - ray.origin[axis]) * ray.inv_direction[axis];

      const KDAccelNode *firstChild = nullptr, *secondChild = nullptr;
      bool belowFirst =
          (ray.origin[axis] < node->SplitPos()) ||
          (ray.origin[axis] == node->SplitPos() && ray.direction[axis] <= 0);
      if (belowFirst) {
        firstChild = node + 1;
        secondChild = &nodes[node->AboveChild()];
      } else {
        firstChild = &nodes[node->AboveChild()];
        secondChild = node + 1;
      }

      if (tPlane > tMax || tPlane <= 0)
        node = firstChild;
      else if (tPlane < tMin)
        node = secondChild;
      else {
        todo[todoPos].node = secondChild;
        todo[todoPos].tMin = tPlane;
        todo[todoPos].tMax = tMax;
        ++todoPos;
        node = firstChild;
        tMax = tPlane;
      }

    } else {
      int curTotalObjects = node->TotalObjects();
      if (curTotalObjects == 1) {
        const ObjectPtr &ptr = objects[node->objectIndex];
        // Check one primitive inside leaf node
        if (ptr->intersect(ray, interaction)) hit = true;
      } else {
        for (int i = 0; i < curTotalObjects; ++i) {
          int index = objectIndices[node->objectIndicesOffset + i];
          const ObjectPtr &ptr = objects[index];
          // Check one primitive inside leaf node
          if (ptr->intersect(ray, interaction)) hit = true;
        }
      }

      // Grab next node to process from todo list
      if (todoPos > 0) {
        --todoPos;
        node = todo[todoPos].node;
        tMin = todo[todoPos].tMin;
        tMax = todo[todoPos].tMax;
      } else
        break;
    }
  }
  return hit;
}