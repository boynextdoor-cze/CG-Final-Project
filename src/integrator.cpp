#include "integrator.h"
#include "utils.h"
#include <omp.h>

#include <utility>

Integrator::Integrator(std::shared_ptr<Camera> cam,
                       std::shared_ptr<Scene> scene, int spp, int max_depth)
    : camera(std::move(cam)), scene(std::move(scene)), spp(spp), max_depth(max_depth) {
}

void Integrator::render() const {
	Vec2i resolution = camera->getImage()->getResolution();
	int cnt = 0;
	Sampler sampler;
// clang-format off
	#pragma omp parallel for schedule(dynamic), shared(cnt, resolution), private(sampler), default(none)
	// clang-format on
	for (int dx = 0; dx < resolution.x(); dx++) {
// clang-format off
		#pragma omp atomic
		// clang-format on
		++cnt;
		printf("\r%.02f%%", cnt * 100.0 / resolution.x());
		sampler.setSeed(omp_get_thread_num());
		for (int dy = 0; dy < resolution.y(); dy++) {
			Vec3f L(0, 0, 0);
			std::vector<Vec2f> ray_samples = camera->getRaySamples((float) dx, (float) dy, spp);
			for (auto &ray_sample: ray_samples) {
				Ray ray = camera->generateRay(ray_sample.x(), ray_sample.y());
				L += radiance(ray, sampler, 0) / (float) ray_samples.size();
			}
			camera->getImage()->setPixel(dx, dy, L);
		}
	}
}

bool isInvalid(Vec3f v) {
	bool isnan = std::isnan(v.x()) || std::isnan(v.y()) || std::isnan(v.z());
	bool isinf = std::isinf(v.x()) || std::isinf(v.y()) || std::isinf(v.z());
	bool toolarge = v.x() > 1.f || v.y() > 1.f || v.z() > 1.f;
	return isnan || isinf || toolarge;
}

Vec3f Integrator::radiance(Ray &ray, Sampler &sampler, int depth) const {
	if (depth >= max_depth) {
		return Vec3f::Zero();
	}
	Interaction in;
	bool scene_intersect = scene->intersect(ray, in);
	if (!scene_intersect || in.type == Interaction::Type::NONE) {
		return Vec3f::Zero();
	}
	if (in.type == Interaction::Type::LIGHT) {
		std::shared_ptr<Light> light = scene->getLight();
		Vec3f L = light->emission(in.normal, ray.direction);
		return L;
	}
	if (in.type == Interaction::Type::GEOMETRY) {
		Vec3f directLight = directLighting(in, sampler);
		directLight = isInvalid(directLight) ? Vec3f::Zero() : directLight;
		Vec3f indirectLight(0, 0, 0);
		in.wo = -ray.direction;
		if (!in.material->isDelta()) {
			Vec3f wi = in.material->sample(in, sampler);
			Ray nextRay(in.pos, wi);
			if (!scene->isShadowed(nextRay)) {
				return directLight;
			}
			Vec3f L = radiance(nextRay, sampler, depth + 1);
			indirectLight = in.material->evaluate(in).cwiseProduct(L) * wi.dot(in.normal) / in.material->pdf(in);
		} else {
			Vec3f wi = -in.wo + 2 * (in.wo.dot(in.normal)) * in.normal;
			Ray nextRay(in.pos, wi);
			if (!scene->isShadowed(nextRay)) {
				return directLight;
			}
			indirectLight = radiance(nextRay, sampler, depth + 1);
		}
		indirectLight = isInvalid(indirectLight) ? Vec3f::Zero() : indirectLight;
		return directLight + indirectLight;
	}
	return Vec3f::Zero();
}

Vec3f Integrator::directLighting(Interaction &interaction, Sampler &sampler) const {
	Vec3f L(0, 0, 0);
	std::shared_ptr<Light> light = scene->getLight();
	Vec3f sample_pos = light->sample(interaction, nullptr, sampler);
	Vec3f ray_dir = sample_pos - interaction.pos;
	Ray shadowRay(interaction.pos, ray_dir);
	if (!scene->isShadowed(shadowRay)) {
		float cos_theta_i = interaction.normal.dot(ray_dir.normalized());
		float cos_theta_o = light->getNormal().dot(-ray_dir.normalized());
		float dist = (float) std::pow((sample_pos - interaction.pos).norm(), 2);
		Vec3f brdf = interaction.material->evaluate(interaction);
		Vec3f light_radiance = light->emission(sample_pos, ray_dir);
		L = cos_theta_i * cos_theta_o * light_radiance.cwiseProduct(brdf);
		L /= dist;
		L /= light->pdf(interaction, sample_pos);
	}
	return L;
}
