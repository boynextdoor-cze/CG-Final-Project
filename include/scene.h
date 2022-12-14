#ifndef SCENE_H_
#define SCENE_H_

#include <vector>

#include "camera.h"
#include "config.h"
#include "geometry.h"
#include "image.h"
#include "interaction.h"
#include "light.h"

class Scene {
public:
	Scene() = default;
	void addObject(std::shared_ptr<TriangleMesh> &geometry);
	[[nodiscard]] const std::shared_ptr<Light> &getLight() const;
	void setLight(const std::shared_ptr<Light> &new_light);
	bool isShadowed(Ray &shadow_ray);
	bool intersect(Ray &ray, Interaction &interaction);
	void buildBVH();

private:
	std::vector<std::shared_ptr<TriangleMesh>> objects;
	std::vector<std::shared_ptr<Triangle>> triangles;
	std::vector<LinearBVHNode> BVHNodes;
	std::shared_ptr<Light> light;
	BVHNode *bvh = nullptr;
	bool bvhHit(Ray &ray, Interaction &interaction, BVHNode *node);
	bool LinearBVHHit(Ray &ray, Interaction &interaction);
	int findSplit(int start, int end);
	void MortonCode();
	BVHNode *generateHierarchy(int start, int end);
	void LinearizeBVH(BVHNode *node);
	BVHNode *LeafNode(int start, int end);
	static BVHNode *InnerNode(BVHNode *left, BVHNode *right);
};

void initSceneFromConfig(const Config &config, std::shared_ptr<Scene> &scene);

#endif//SCENE_H_
