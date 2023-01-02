#include "light.h"

#include "utils.h"
#include <utility>

constexpr uint32_t SAMPLE_NUM = 16;

Light::Light(Vec3f pos, Vec3f color)
    : position(std::move(pos)), radiance(std::move(color)) {}

SquareAreaLight::SquareAreaLight(const Vec3f &pos, const Vec3f &color,
                                 const Vec2f &size)
    : Light(pos, color), size(size) {
	Vec3f v1, v2, v3, v4;
	v1 = pos + Vec3f(size.x() / 2, 0.f, -size.y() / 2);
	v2 = pos + Vec3f(-size.x() / 2, 0.f, -size.y() / 2);
	v3 = pos + Vec3f(-size.x() / 2, 0.f, size.y() / 2);
	v4 = pos + Vec3f(size.x() / 2, 0.f, size.y() / 2);
	light_mesh = TriangleMesh({v1, v2, v3, v4}, {Vec3f(0, -1, 0)},
	                          {0, 1, 2, 0, 2, 3}, {0, 0, 0, 0, 0, 0});
}

Vec3f SquareAreaLight::emission(const Vec3f &pos, const Vec3f &dir) const {
	return radiance * std::max(dir.dot(Vec3f(0, -1, 0)), 0.f);
}

float SquareAreaLight::pdf(const Interaction &interaction, Vec3f pos) const {
	return 1 / (size.x() * size.y());
}

Vec3f SquareAreaLight::sample(Interaction &interaction,
                              Sampler &sampler) const {
	Vec2f offset = (sampler.get2D() - Vec2f(0.5f, 0.5f)).cwiseProduct(size);
	Vec3f result = position + Vec3f(offset.x(), 0.f, offset.y());
	return result;
}

bool SquareAreaLight::intersect(Ray &ray, Interaction &interaction) const {
	if (light_mesh.intersect(ray, interaction)) {
		interaction.type = Interaction::Type::LIGHT;
		return true;
	}
	return false;
}
