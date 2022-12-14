#include "bsdf.h"
#include "utils.h"

#include <utility>

IdealDiffusion::IdealDiffusion(Vec3f color) : color(std::move(color)) {}

Vec3f IdealDiffusion::evaluate(Interaction &interaction) const {
	return color / PI;
}

float IdealDiffusion::pdf(Interaction &interaction) const {
	float cos_theta = interaction.normal.dot(interaction.wi);
	return cos_theta / PI;
}

Vec3f IdealDiffusion::sample(Interaction &interaction, Sampler &sampler) const {
	Vec2f sample = sampler.get2D();
	float phi = 2.f * PI * sample.y();
	//float theta = std::acos(std::sqrt(1.f - (float) std::pow(sample.x(), 2)));
	float theta = std::acos(std::sqrt(1.f - sample.x()));
	Vec3f hemi_wi = Vec3f(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
	Vec3f hemi_normal(0, 0, 1);
	Mat3f transform = Eigen::Quaternionf::FromTwoVectors(hemi_normal, interaction.normal).toRotationMatrix();
	Vec3f wi = (transform * hemi_wi).normalized();
	interaction.wi = wi;
	return wi;
}
// return whether the bsdf is perfect transparent or perfect reflection
bool IdealDiffusion::isDelta() const {
	return false;
}
IdealReflection::IdealReflection(Vec3f color) {
	this->color = std::move(color);
}
Vec3f IdealReflection::evaluate(Interaction &interaction) const {
	return Vec3f::Zero();
}
float IdealReflection::pdf(Interaction &interaction) const {
	return 0;
}
Vec3f IdealReflection::sample(Interaction &interaction, Sampler &sampler) const {
	return Vec3f::Zero();
}
bool IdealReflection::isDelta() const {
	return true;
}
