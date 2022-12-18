#include "bsdf.h"
#include "utils.h"

#include <utility>

IdealDiffusion::IdealDiffusion(const Vec3f &color) : color(color) {}

Vec3f IdealDiffusion::evaluate(Interaction &interaction) const {
  // why 1 / pi ?
  constexpr float DIFFUSE_BRDF = INV_PI;
  return color * DIFFUSE_BRDF;
}

float IdealDiffusion::pdf(Interaction &interaction) const {
  // All direction in interaction should be normalized
  float cosine_theta = interaction.wi.dot(interaction.normal);
  return cosine_theta * INV_PI;
}

float IdealDiffusion::sample(Interaction &interaction, Sampler &sampler) const {
  float ksi_1 = sampler.get1D(), ksi_2 = sampler.get1D();
  // cosine weighted sample
  float x = std::cos(2 * PI * ksi_2) * std::sqrt(ksi_1);
  float y = std::sin(2 * PI * ksi_2) * std::sqrt(ksi_1);
  float z = std::sqrt(1 - ksi_1);
  Vec3f local_wi = Vec3f(x, y, z);
  Mat3f rotate = Mat3f::Zero();
  rotate =
      Eigen::Quaternionf::FromTwoVectors(Vec3f(0, 0, 1.0f), interaction.normal)
          .toRotationMatrix();
  interaction.wi = (rotate * local_wi).normalized();
  return pdf(interaction);
}
/// return whether the bsdf is perfect transparent or perfect reflection
bool IdealDiffusion::isDelta() const { return false; }

IdealSpecular::IdealSpecular(const Vec3f &color) : color(color) {}

Vec3f IdealSpecular::evaluate(Interaction &interaction) const {
  Vec3f ideal_wi =
      -interaction.wo +
      2 * interaction.normal.dot(interaction.wo) * interaction.normal;
  if ((ideal_wi - interaction.wi).dot(ideal_wi - interaction.wi) < EPS * EPS) {
    return color / std::fabs(interaction.normal.dot(interaction.wi));
  } else
    return Vec3f(0, 0, 0);
}

float IdealSpecular::pdf(Interaction &interaction) const { return 1.0f; }

float IdealSpecular::sample(Interaction &interaction, Sampler &sampler) const {
  interaction.wi =
      -interaction.wo +
      2 * interaction.normal.dot(interaction.wo) * interaction.normal;
  return pdf(interaction);
}

bool IdealSpecular::isDelta() const { return true; }