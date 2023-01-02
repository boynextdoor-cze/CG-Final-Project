#include <iostream>
#include <chrono>

#include "integrator.h"
#include "config_io.h"
#include "config.h"

#include <fstream>

int main(int argc, char *argv[]) {
  /// load config from json file
  Config config;
  std::ifstream fin;
  std::string file_name;
  if (argc == 1) {
    std::cout << "No json specified, use default path." << std::endl;
    file_name = "configs/simple.json";
    fin.open("configs/simple.json");
  } else {
    file_name = argv[1];
    fin.open(argv[1]);
  }
  if (!fin.is_open()) {
    std::cerr << "Can not open json file. Exit." << std::endl;
    exit(0);
  } else {
    std::cout << "Json file loaded from " << file_name << std::endl;
  }
  // parse json object to Config
  try {
    nlohmann::json j;
    fin >> j;
    nlohmann::from_json(j, config);
    fin.close();
  } catch (nlohmann::json::exception &ex) {
    fin.close();
    std::cerr << "Error:" << ex.what() << std::endl;
    exit(-1);
  }
  std::cout << "Parsed json to config. Start building scene..." << std::endl;
  // initialize all settings from config
  // set image resolution.
  std::shared_ptr<ImageRGB> rendered_img
      = std::make_shared<ImageRGB>(config.image_resolution[0], config.image_resolution[1]);
  std::cout << "Image resolution: "
            << config.image_resolution[0] << " x " << config.image_resolution[1] << std::endl;
  // set camera
  std::shared_ptr<Camera> camera = std::make_shared<Camera>(config.cam_config, rendered_img);
  // construct scene.
  auto scene = std::make_shared<Scene>();
  initSceneFromConfig(config, scene);
  // init integrator
  std::unique_ptr<Integrator> integrator
      = std::make_unique<Integrator>(camera, scene, config.spp, config.max_depth);
  std::cout << "Start Rendering..." << std::endl;
  auto start = std::chrono::steady_clock::now();
  // render scene
  integrator->render();
  auto end = std::chrono::steady_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  std::cout << "\nRender Finished in " << time << "s." << std::endl;
  rendered_img->writeImgToFile("result.png");
  std::cout << "Image saved to disk." << std::endl;
  return 0;
}