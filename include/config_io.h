#ifndef CONFIG_IO_H_
#define CONFIG_IO_H_

#define JSON_USE_IMPLICIT_CONVERSIONS 0

#include <nlohmann/json.hpp>

#include "config.h"

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config::LightConfig, position, size,
                                   radiance);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config::CamConfig, position, look_at, ref_up,
                                   vertical_fov, focal_length);

// add your own bsdf name if needed
NLOHMANN_JSON_SERIALIZE_ENUM(MaterialType,
                             {{MaterialType::DIFFUSE, "diffuse"},
                              {MaterialType::SPECULAR, "specular"}});

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config::MaterialConfig, color, type, name);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config::ObjConfig, obj_file_path,
                                   material_name, translate, scale, has_bvh);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config::NURBSConfig, nurbs_json_path,
                                   material_name, translate, scale);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Config, spp, max_depth, image_resolution,
                                   cam_config, light_config, materials, objects,
                                   nurbs);

// NURBS Coniguration
NLOHMANN_JSON_SERIALIZE_ENUM(NURBSShapeType,
                             {{NURBSShapeType::SURFACE, "surface"}});

NLOHMANN_JSON_SERIALIZE_ENUM(TrimType, {{TrimType::SPLINE, "spline"},
                                        {TrimType::CONTAINER, "container"}});

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig::ControlPoint, points, weights);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig::TrimCurveData, control_points,
                                   degree, dimension, knotvector, rational,
                                   reversed, type);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig::InnerTrimData, count, data,
                                   reversed, type);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig::TrimData, count, data);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig::NURBSData, control_points,
                                   degree_u, degree_v, dimension, knotvector_u,
                                   knotvector_v, rational, reversed, size_u,
                                   size_v, trims);

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(NURBSConfig, count, data, type);
#endif// CONFIG_IO_H_