#include "scene.h"
#include "load_obj.h"
#include "bezier.h"

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

  {
    /* Test: add NURBSes into scene. */
    std::shared_ptr<BSDF> GREEN, YELLOW;
    GREEN = std::make_shared<IdealDiffusion>(Vec3f(0.04, 0.7, 0.04));
    YELLOW = std::make_shared<IdealDiffusion>(Vec3f(0.7, 0.7, 0.04));


    // set duck_1
    std::shared_ptr<NURBS> tmpDuck1 = std::make_shared<NURBS>(13, 12, 4, 4);
    tmpDuck1->setKnotM({-1.5708, -1.5708, -1.5708, -1.5708, -1.0472, -0.523599, 0, 0.523599, 0.808217,
                       1.04015, 1.0472, 1.24824, 1.29714, 1.46148, 1.5708, 1.5708, 1.5708, 1.5708});
    tmpDuck1->setKnotN({-3.14159, -3.14159, -3.14159, -3.14159, -2.61799, -2.0944, -1.0472, -0.523599,
                       6.66134e-016, 0.523599, 1.0472, 2.0944, 2.61799, 3.14159, 3.14159, 3.14159, 3.14159});
    tmpDuck1->setControlPoint(readControlPoints(R"(../assets/duck1.ctrlpts)", 13, 12));
    tmpDuck1->setWeight(readWeights(R"(../assets/duck1.weights)", 13, 12));
		tmpDuck1->refine();
		std::shared_ptr<TriangleMesh> tmpObject1 = tmpDuck1->generateMesh(Uniform);
		tmpObject1->setMaterial(GREEN);
    tmpObject1->buildBVH();

    // set duck_2
    std::shared_ptr<NURBS> tmpDuck2 = std::make_shared<NURBS>(8, 9, 4, 4);
    tmpDuck2->setKnotM({0, 0, 0, 0, 0.145456, 0.265731, 0.436096, 0.583258, 0.847704, 1, 1, 1, 1});
    tmpDuck2->setKnotN({0, 0, 0, 0, 0.179541, 0.317924, 0.485586, 0.507528, 0.709398, 0.813231, 1, 1, 1, 1});
    tmpDuck2->setControlPoint(readControlPoints(R"(../assets/duck2.ctrlpts)", 8, 9));
		tmpDuck2->refine();
		std::shared_ptr<TriangleMesh> tmpObject2 = tmpDuck2->generateMesh(Uniform);
    tmpObject2->setMaterial(YELLOW);
    tmpObject2->buildBVH();

    // set duck_3
    std::shared_ptr<NURBS> tmpDuck3 = std::make_shared<NURBS>(5, 5, 4, 4);
    tmpDuck3->setKnotM({0, 0, 0, 0, 0.333333, 0.666667, 1, 1, 1, 1});
    tmpDuck3->setKnotN({0, 0, 0, 0, 0.333333, 0.666667, 1, 1, 1, 1});
    tmpDuck3->setControlPoint(readControlPoints(R"(../assets/duck3.ctrlpts)", 5, 5));
		tmpDuck3->refine();
    std::shared_ptr<TriangleMesh> tmpObject3 = tmpDuck3->generateMesh(Uniform);
    tmpObject3->setMaterial(YELLOW);
    tmpObject3->buildBVH();		

    // add objects to scene
    scene->addObject(tmpObject1);
    scene->addObject(tmpObject2);
    scene->addObject(tmpObject3);
  }

  /* Build Global BVH. */
  scene->buildBVH();
}