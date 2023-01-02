#ifndef SCENE_H_
#define SCENE_H_

#include <vector>

#include "accel.h"
#include "camera.h"
#include "config.h"
#include "geometry.h"
#include "image.h"
#include "interaction.h"
#include "light.h"

class Scene {
public:
	Scene() = default;
	void addObject(const ObjectPtr &geometry);
	[[nodiscard]] const std::shared_ptr<Light> &getLight() const;
	void setLight(const std::shared_ptr<Light> &new_light);
	bool isShadowed(Ray &shadow_ray, float light_dist = (float) 1e10);
	bool intersect(Ray &ray, Interaction &interaction);
	void buildBVH();

private:
	std::vector<ObjectPtr> objects;
	std::shared_ptr<Light> light;
	BVHAccelPtr bvh;
};

void initSceneFromConfig(const Config &config, std::shared_ptr<Scene> &scene);

#endif// SCENE_H_
