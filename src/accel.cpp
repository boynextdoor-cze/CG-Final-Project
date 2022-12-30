#include "accel.h"
#include "utils.h"
#include <iostream>
#include <stdlib.h>
#include <ctime>

#undef NDEBUG
#include <cassert>

struct BVHObjectInfo {
  int index;
  Bounds3 bound;
  Vec3f centroid;
  BVHObjectInfo(int _index, const Bounds3 &_bound)
      : index(_index), bound(_bound), centroid(_bound.Centroid()) {}
  BVHObjectInfo() = default;
};

struct BucketInfo {
  int count = 0;
  Bounds3 bounds;
};

BVHAccel::~BVHAccel() {
  if (nodes != nullptr)
    free(nodes);
}

BVHAccel::BVHAccel(std::vector<ObjectPtr> &p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      objects(std::move(p)) {
  time_t start, stop;
  time(&start);

  if (objects.empty())
    return;

  int total_nodes = 0;

  std::vector<BVHObjectInfo> object_info(objects.size());
  std::vector<ObjectPtr> ordered_objects;
  ordered_objects.reserve(objects.size());

  for (int i = 0; i < objects.size(); ++i)
    object_info[i] = {i, objects[i]->getBounds()};

  root = recursiveBuild(object_info, 0, objects.size(), total_nodes,
                        ordered_objects);

  recursiveCheck(root);

  objects.swap(ordered_objects);
  object_info.resize(0);

  // int treeBytes = total_nodes * sizeof(LinearBVHNode) + sizeof(*this) +
  //                 objects.size() * sizeof(objects[0]);

  // Cache aligned, assume cache line is 64 bytes
//  nodes =
//      (LinearBVHNode *)aligned_alloc(64, total_nodes * sizeof(LinearBVHNode));
  nodes = (LinearBVHNode *)malloc(total_nodes * sizeof(LinearBVHNode));
  // nodes = std::vector<LinearBVHNode>(total_nodes);
  int offset = 0;
  flattenBVHTree(root, offset);
  assert(offset == total_nodes);
  // assert (false);
  for (int i = 0; i < total_nodes; i++) {
    auto &node = nodes[i];
    if (node.objects_num > 0) {
      assert(node.object_offset < objects.size() && node.object_offset >= 0);
      // assert(node.objects_num <= 2);
      assert(node.objects_num + node.object_offset <= objects.size());
    } else {
      assert(node.right_offset >= 0 && node.right_offset < total_nodes);
      assert(node.split_axis >= 0 && node.split_axis < 3);
    }
  }

  if (offset != total_nodes) {
    std::cout << "Error" << std::endl;
    exit(1);
  }
  time(&stop);
  double diff = difftime(stop, start);
  int hrs = (int)diff / 3600;
  int mins = ((int)diff / 60) - (hrs * 60);
  int secs = (int)diff - (hrs * 3600) - (mins * 60);

  printf(
      "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
      hrs, mins, secs);
}

BVHNodePtr BVHAccel::recursiveBuild(std::vector<BVHObjectInfo> &object_info,
                                    int begin, int end, int &total_nodes,
                                    std::vector<ObjectPtr> &ordered_obj) {
  if (begin >= end) {
    std::cout << "object empty" << std::endl;
    return nullptr;
  }

  BVHNodePtr node = std::make_shared<BVHNode>();
  total_nodes++;
  int objects_num = end - begin;
  Bounds3 bounds;
  for (int i = begin; i < end; ++i)
    bounds = Union(bounds, object_info[i].bound);
  // std::cout << "bounds' area is " << bounds.SurfaceArea() << std::endl;

  if (objects_num == 1) {
    int first = ordered_obj.size();
    for (int i = begin; i < end; i++) {
      int index = object_info[i].index;
      ordered_obj.push_back(objects[index]);
    }
    node->initLeaf(first, objects_num, bounds);
    return node;
  } else {
    Bounds3 centroidBounds;
    for (int i = begin; i < end; ++i)
      centroidBounds = Union(centroidBounds, object_info[i].centroid);
    int dim = centroidBounds.maxExtent();
    int mid = (begin + end) / 2;

    if (std::fabs(centroidBounds.pMax[dim] - centroidBounds.pMin[dim]) < EPS * EPS) {
      int first = ordered_obj.size();
      for (int i = begin; i < end; i++) {
        int index = object_info[i].index;
        ordered_obj.push_back(objects[index]);
      }
      node->initLeaf(first, objects_num, bounds);
      return node;
    } else {
      // Only SAH
      if (objects_num <= 2) {
        mid = (begin + end) / 2;
        std::nth_element(&object_info[begin], &object_info[mid],
                         &object_info[end - 1] + 1,
                         [dim](const BVHObjectInfo &a, const BVHObjectInfo &b) {
                           return a.centroid[dim] < b.centroid[dim];
                         });
      } else {
        const int B = 32;
        BucketInfo buckets[B];
        // determine which bucket each object is in
        for (int i = begin; i < end; i++) {
          int b = B * centroidBounds.Offset(object_info[i].centroid)[dim];
          // on the border
          if (b == B)
            b = B - 1;
          if (b < 0) {
            std::cout << "Something wrong" << std::endl;
          }
          buckets[b].count++;
          buckets[b].bounds = Union(buckets[b].bounds, object_info[i].bound);
        }
        float cost[B - 1];

        Bounds3 pre, suf;
        int pre_count = 0, suf_count = 0;
        for (int i = 0; i < B - 1; i++) {
          pre = Union(pre, buckets[i].bounds);
          pre_count += buckets[i].count;
          cost[i] = 1.f + pre_count * pre.SurfaceArea() / bounds.SurfaceArea();
        }
        for (int i = B - 1; i >= 1; i--) {
          suf = Union(suf, buckets[i].bounds);
          suf_count += buckets[i].count;
          cost[i - 1] += suf_count * suf.SurfaceArea() / bounds.SurfaceArea();
        }

        float min_cost = cost[0];
        int min_split_index = 0;
        for (int i = 1; i < B - 1; i++) {
          if (cost[i] < min_cost) {
            min_cost = cost[i];
            min_split_index = i;
          }
        }

        float leaf_cost = objects_num;
        if (objects_num > maxPrimsInNode || min_cost < leaf_cost) {
          // should not make this a leaf node
          BVHObjectInfo *pmid =
              std::partition(&object_info[begin], &object_info[end - 1] + 1,
                             [=](const BVHObjectInfo &pi) {
                               int b =
                                   B * centroidBounds.Offset(pi.centroid)[dim];
                               if (b == B)
                                 b = B - 1;
                               return b <= min_split_index;
                             });
          mid = pmid - &object_info[0];

        } else {
          int first = ordered_obj.size();
          for (int i = begin; i < end; i++) {
            int index = object_info[i].index;
            ordered_obj.push_back(objects[index]);
          }
          node->initLeaf(first, objects_num, bounds);
          return node;
        }
      }
      node->initInterior(
          dim,
          recursiveBuild(object_info, begin, mid, total_nodes, ordered_obj),
          recursiveBuild(object_info, mid, end, total_nodes, ordered_obj));
    }
  }
  return node;
}

void BVHAccel::recursiveCheck(BVHNodePtr node) {
  if (node->objects_num > 0) {
    if (node->left != nullptr || node->right != nullptr) {
      std::cout << "Wrong, leaf node cannot have child." << std::endl;
    }
    if (!(node->objects_num + node->first_object_offset <= objects.size())) {
      std::cout << "offset is " << node->first_object_offset
                << ", total object is " << node->objects_num << std::endl;
    }
    assert(node->objects_num + node->first_object_offset <= objects.size());
    return;
  } else {
    if (node->left == nullptr || node->right == nullptr ||
        node->objects_num != 0) {
      std::cout << "Wrong, interior node must have two children and "
                   "objects_num equal to zero."
                << std::endl;
      std::cout << node->left << " " << node->right << " " << node->objects_num
                << std::endl;
    }
    recursiveCheck(node->left);
    recursiveCheck(node->right);
  }
}

int BVHAccel::flattenBVHTree(BVHNodePtr node, int &offset) {
  LinearBVHNode *linearNode = (nodes + offset);
  linearNode->bounds = node->bounds;
  int myOffset = offset;
  offset++;
  if (node->objects_num > 0) {
    if (node->objects_num > 100) {
      std::cout << node->objects_num << std::endl;
      exit(1);
    }
    linearNode->object_offset = node->first_object_offset;
    linearNode->objects_num = node->objects_num;
    assert(node->first_object_offset + node->objects_num <= objects.size());
    assert(linearNode->object_offset + linearNode->objects_num <=
           objects.size());
  } else {
    // Create interior flattened BVH node
    linearNode->split_axis = node->split_axis;
    linearNode->objects_num = 0;
    flattenBVHTree(node->left, offset);
    linearNode->right_offset = flattenBVHTree(node->right, offset);
    // assert(linearNode->right_offset < nodes.size());
  }
  return myOffset;
}

bool BVHAccel::getIntersection(const Ray &ray, Interaction &optimal,
                               bool hitAny) const {
  if (!nodes)
    return false;
  bool hit = false;
  Vec3f inv_direction(1 / ray.direction.x(), 1 / ray.direction.y(),
                      1 / ray.direction.z());
  int dirIsNeg[3] = {inv_direction.x() < 0, inv_direction.y() < 0,
                     inv_direction.z() < 0};
  // Follow ray through BVH nodes to find primitive intersections
  int toVisitOffset = 0, currentNodeIndex = 0;
  int nodesToVisit[64];
  while (true) {
    const LinearBVHNode *node = &nodes[currentNodeIndex];
    auto [flag, t_enter] = node->bounds.IntersectWithTime(ray);
    if (flag && t_enter < optimal.dist) {
      if (node->objects_num > 0) {
        for (int i = 0; i < node->objects_num; ++i)
          if (objects[node->object_offset + i]->intersect(ray, optimal)) {
            hit = true;
            if (hitAny) return true;
          }
        if (toVisitOffset == 0)
          break;
        currentNodeIndex = nodesToVisit[--toVisitOffset];
      } else {
        // Put far BVH node on _nodesToVisit_ stack, advance to near
        // node
        if (dirIsNeg[node->split_axis]) {
          nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
          currentNodeIndex = node->right_offset;
        } else {
          nodesToVisit[toVisitOffset++] = node->right_offset;
          currentNodeIndex = currentNodeIndex + 1;
        }
      }
    } else {
      if (toVisitOffset == 0)
        break;
      currentNodeIndex = nodesToVisit[--toVisitOffset];
    }
  }
  return hit;
}

void BVHAccel::getIntersection(BVHNodePtr node, const Ray &ray,
                               Interaction &optimal) const {
  if (!node)
    return;
  auto [flag, t_enter] = node->bounds.IntersectWithTime(ray);
  if (!(flag && t_enter < optimal.dist))
    return;
  if (node->objects_num > 0) {
    Interaction inter;
    for (int i = 0; i < node->objects_num; i++) {
      objects[node->first_object_offset + i]->intersect(ray, inter);
    }
    if (inter.dist < optimal.dist)
      optimal = inter;
    return;
  }
  BVHNodePtr l = node->left, r = node->right;
  auto [flag_l, t_enter_l] = l->bounds.IntersectWithTime(ray);
  auto [flag_r, t_enter_r] = r->bounds.IntersectWithTime(ray);
  if (t_enter_l > t_enter_r) {
    std::swap(flag_l, flag_r);
    std::swap(t_enter_l, t_enter_r);
    std::swap(l, r);
  }
  if (flag_l && t_enter_l < optimal.dist)
    getIntersection(l, ray, optimal);
  if (flag_r && t_enter_r < optimal.dist)
    getIntersection(r, ray, optimal);
}

Bounds3 BVHAccel::WorldBound() const { return root->bounds; }

BVHNode::BVHNode() {
  bounds = Bounds3();
  left = nullptr;
  right = nullptr;
  split_axis = 0;
  first_object_offset = 0;
  objects_num = 0;
}

void BVHNode::initInterior(int _split_axis, BVHNodePtr _left,
                           BVHNodePtr _right) {
  left = _left;
  right = _right;
  split_axis = _split_axis;
  bounds = Union(left->bounds, right->bounds);
  objects_num = 0;
  first_object_offset = 0;
}

void BVHNode::initLeaf(int first, int n, const Bounds3 &bound) {
  objects_num = n;
  first_object_offset = first;
  bounds = bound;
  left = nullptr;
  right = nullptr;
  split_axis = 0;
}