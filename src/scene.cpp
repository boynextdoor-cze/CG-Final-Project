#include "scene.h"
#include "accel.h"
#include "geometry.h"
#include "load_obj.h"

#include <iostream>
#include <stack>
#include <utility>

void Scene::addObject(std::shared_ptr<TriangleMesh> &mesh) {
	objects.push_back(mesh);
}

void Scene::setLight(const std::shared_ptr<Light> &new_light) {
	light = new_light;
}
bool Scene::isShadowed(Ray &shadow_ray) {
	Interaction in;
	return intersect(shadow_ray, in) && in.type == Interaction::Type::GEOMETRY;
}

bool Scene::intersect(Ray &ray, Interaction &interaction) {
	if (bvh != nullptr) {
		light->intersect(ray, interaction);
		//bool hit = bvhHit(ray, interaction, bvh);
		bool hit = LinearBVHHit(ray, interaction);
		return hit;
	} else {
		light->intersect(ray, interaction);
		for (const auto &obj: objects) {
			Interaction cur_it;
			if (obj->intersect(ray, cur_it) && (cur_it.dist < interaction.dist)) {
				interaction = cur_it;
			}
		}
		return interaction.type != Interaction::Type::NONE;
	}
}

const std::shared_ptr<Light> &Scene::getLight() const {
	return light;
}
bool cmp(const std::shared_ptr<Triangle> &a, const std::shared_ptr<Triangle> &b) {
	return a->morton_code < b->morton_code;
}
void Scene::MortonCode() {
	AABB scene_box = objects[0]->generateAABB();
	for (const auto &obj: objects) {
		AABB obj_aabb = obj->generateAABB();
		scene_box = AABB(scene_box, obj_aabb);
	}
	for (const auto &obj: objects) {
		auto tris = obj->generateMortonCode(scene_box);
		triangles.insert(triangles.end(), tris.begin(), tris.end());
	}
	std::sort(triangles.begin(), triangles.end(), cmp);
}
void Scene::buildBVH() {
	MortonCode();
	bvh = generateHierarchy(0, (int) triangles.size() - 1);
	LinearizeBVH(bvh);
}

BVHNode *Scene::generateHierarchy(int start, int end) {
	int len = end - start + 1;
	if (len >= 1 && len <= 8) {
		return LeafNode(start, end);
	} else {
		int split = findSplit(start, end);
		BVHNode *left = generateHierarchy(start, split);
		BVHNode *right = generateHierarchy(split + 1, end);
		return InnerNode(left, right);
	}
}

BVHNode *Scene::LeafNode(int start, int end) {
	auto *node = new BVHNode();
	node->start = start;
	node->end = end;
	AABB bbox;
	for (int i = start; i <= end; i++) {
		AABB tri_bbox = triangles[i]->aabb;
		if (i == start) {
			bbox = tri_bbox;
		} else {
			bbox = AABB(bbox, tri_bbox);
		}
	}
	node->aabb = bbox;
	node->left = nullptr;
	node->right = nullptr;
	node->size = 1;
	return node;
}

BVHNode *Scene::InnerNode(BVHNode *left, BVHNode *right) {
	auto *node = new BVHNode();
	node->aabb = AABB(left->aabb, right->aabb);
	node->left = left;
	node->right = right;
	node->size = left->size + right->size + 1;
	return node;
}

void Scene::LinearizeBVH(BVHNode *node) {
	if (node == nullptr) {
		return;
	}
	LinearBVHNode linearNode;
	linearNode.aabb = node->aabb;
	if (node->left == nullptr && node->right == nullptr) {
		linearNode.start = node->start;
		linearNode.end = node->end;
	} else {
		linearNode.right = (int) BVHNodes.size() + node->left->size + 1;
	}
	BVHNodes.push_back(linearNode);
	if (node->left != nullptr) {
		LinearizeBVH(node->left);
	}
	if (node->right != nullptr) {
		LinearizeBVH(node->right);
	}
}

bool Scene::bvhHit(Ray &ray, Interaction &interaction, BVHNode *node) {
	float t_in, t_out;
	if (node == nullptr || !node->aabb.intersect(ray, t_in, t_out)) {
		return false;
	}
	if (node->left == nullptr && node->right == nullptr) {
		for (int i = node->start; i <= node->end; i++) {
			Interaction cur_inter;
			if (triangles[i]->intersect(ray, cur_inter) && cur_inter.dist < interaction.dist) {
				interaction = cur_inter;
			}
		}
		return interaction.type != Interaction::Type::NONE;
	} else {
		bool hit_left = false;
		bool hit_right = false;
		if (node->left != nullptr)
			hit_left = bvhHit(ray, interaction, node->left);
		if (node->right != nullptr)
			hit_right = bvhHit(ray, interaction, node->right);
		return hit_left || hit_right;
	}
}

bool Scene::LinearBVHHit(Ray &ray, Interaction &interaction) {
	std::stack<int> nodesToVisit;
	int cur_node = 0;
	bool hit = false;
	while (true) {
		LinearBVHNode node = BVHNodes[cur_node];
		float t_in, t_out;
		if (!node.aabb.intersect(ray, t_in, t_out)) {
			if (nodesToVisit.empty()) {
				break;
			}
			cur_node = nodesToVisit.top();
			nodesToVisit.pop();
		}
		if (node.start != -1) {
			for (int i = node.start; i <= node.end; i++) {
				Interaction cur_inter;
				if (triangles[i]->intersect(ray, cur_inter) && cur_inter.dist < interaction.dist) {
					interaction = cur_inter;
				}
			}
			if (interaction.type != Interaction::Type::NONE)
				hit = true;
			if (nodesToVisit.empty()) {
				break;
			}
			cur_node = nodesToVisit.top();
			nodesToVisit.pop();
		} else {
			float t_in_left, t_out_left;
			float t_in_right, t_out_right;
			bool hit_right = BVHNodes[node.right].aabb.intersect(ray, t_in_right, t_out_right);
			bool hit_left = BVHNodes[cur_node + 1].aabb.intersect(ray, t_in_left, t_out_left);
			if (hit_left && hit_right) {
				if (t_in_left < t_in_right) {
					nodesToVisit.push(node.right);
					cur_node++;
				} else {
					nodesToVisit.push(cur_node + 1);
					cur_node = node.right;
				}
			} else if (hit_left) {
				cur_node++;
			} else if (hit_right) {
				cur_node = node.right;
			} else {
				if (nodesToVisit.empty()) {
					break;
				}
				cur_node = nodesToVisit.top();
				nodesToVisit.pop();
			}
		}
	}
	return hit;
}

int Scene::findSplit(int start, int end) {
	std::shared_ptr<Triangle> firstTriangle = triangles[start];
	std::shared_ptr<Triangle> lastTriangle = triangles[end];
	unsigned int firstCode = firstTriangle->morton_code;
	unsigned int lastCode = lastTriangle->morton_code;
	if (firstCode == lastCode) {
		return (start + end) >> 1;
	}
	int commonPrefix = __builtin_clz(firstCode ^ lastCode);
	int split = start;
	int step = end - start;
	do {
		step = (step + 1) >> 1;
		int newSplit = split + step;
		if (newSplit < end) {
			unsigned int splitCode = triangles[newSplit]->morton_code;
			int splitPrefix = __builtin_clz(firstCode ^ splitCode);
			if (splitPrefix > commonPrefix) {
				split = newSplit;
			}
		}
	} while (step > 1);
	return split;
}

void initSceneFromConfig(const Config &config, std::shared_ptr<Scene> &scene) {
	// add square light to scene.
	std::shared_ptr<Light> light = std::make_shared<SquareAreaLight>(Vec3f(config.light_config.position),
	                                                                 Vec3f(config.light_config.radiance),
	                                                                 Vec2f(config.light_config.size));
	scene->setLight(light);
	// init all materials.
	std::map<std::string, std::shared_ptr<BSDF>> mat_list;
	for (const auto &mat: config.materials) {
		std::shared_ptr<BSDF> p_mat;
		switch (mat.type) {
			case MaterialType::DIFFUSE: {
				p_mat = std::make_shared<IdealDiffusion>(Vec3f(mat.color));
				mat_list[mat.name] = p_mat;
				break;
			}
			case MaterialType::SPECULAR: {
				p_mat = std::make_shared<IdealReflection>(Vec3f(mat.color));
				mat_list[mat.name] = p_mat;
				break;
			}
			default: {
				std::cerr << "unsupported material type!" << std::endl;
				exit(-1);
			}
		}
	}
	// add mesh objects to scene. Translation and scaling are directly applied to vertex coordinates.
	// then set corresponding material by name.
	std::cout << "loading obj files..." << std::endl;
	for (auto &object: config.objects) {
		auto mesh_obj = makeMeshObject(object.obj_file_path, Vec3f(object.translate), object.scale);
		mesh_obj->setMaterial(mat_list[object.material_name]);
		scene->addObject(mesh_obj);
	}
	scene->buildBVH();
}