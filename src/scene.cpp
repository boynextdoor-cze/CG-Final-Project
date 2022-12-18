#include "scene.h"
#include "load_obj.h"

#include <iostream>
#include <utility>

void Scene::addObject(std::shared_ptr<TriangleMesh> &mesh) {
  objects.push_back(mesh);
  // mesh->pushAllTriangles(objects);
}

void Scene::setLight(const std::shared_ptr<Light> &new_light) {
  light = new_light;
}
bool Scene::isShadowed(Ray &shadow_ray, float light_dist) {
  Interaction in;
  if (bvh != nullptr) {
    in.dist = light_dist;
    return bvh->getIntersection(shadow_ray, in, true);
  } else
    return intersect(shadow_ray, in) && in.type == Interaction::Type::GEOMETRY;
}

bool Scene::intersect(Ray &ray, Interaction &interaction) {
  light->intersect(ray, interaction);
  if (bvh != nullptr) {
    // bvh->getIntersection(bvh->root, ray, interaction);
    bvh->getIntersection(ray, interaction);
  } else {
    for (const auto &obj : objects) {
      Interaction cur_it;
      if (obj->intersect(ray, cur_it) && (cur_it.dist < interaction.dist)) {
        interaction = cur_it;
      }
    }
  }
  return interaction.type != Interaction::Type::NONE;
}

const std::shared_ptr<Light> &Scene::getLight() const { return light; }

void Scene::buildBVH() {
  bvh = std::make_shared<BVHAccel>(objects);
}

void initSceneFromConfig(const Config &config, std::shared_ptr<Scene> &scene) {
  // add square light to scene.
  std::shared_ptr<Light> light = std::make_shared<SquareAreaLight>(
      Vec3f(config.light_config.position), Vec3f(config.light_config.radiance),
      Vec2f(config.light_config.size));
  scene->setLight(light);
  // init all materials.
  std::map<std::string, std::shared_ptr<BSDF>> mat_list;
  for (const auto &mat : config.materials) {
    std::shared_ptr<BSDF> p_mat;
    switch (mat.type) {
    case MaterialType::DIFFUSE: {
      p_mat = std::make_shared<IdealDiffusion>(Vec3f(mat.color));
      mat_list[mat.name] = p_mat;
      break;
    }
    case MaterialType::SPECULAR: {
      p_mat = std::make_shared<IdealSpecular>(Vec3f(mat.color));
      mat_list[mat.name] = p_mat;
      break;
    }
    default: {
      std::cerr << "unsupported material type!" << std::endl;
      exit(-1);
    }
    }
  }
  // add mesh objects to scene. Translation and scaling are directly applied to
  // vertex coordinates. then set corresponding material by name.
  std::cout << "loading obj files..." << std::endl;
  for (auto &object : config.objects) {
    auto mesh_obj = makeMeshObject(object.obj_file_path,
                                   Vec3f(object.translate), object.scale);
    mesh_obj->setMaterial(mat_list[object.material_name]);
    if (object.has_bvh) {
      mesh_obj->buildBVH();
    }
    scene->addObject(mesh_obj);
  }

  /* Build Global BVH. */
  scene->buildBVH();
}