#include "integrator.h"
#include "utils.h"
#include <omp.h>

#include <chrono>
#include <iostream>

Integrator::Integrator(std::shared_ptr<Camera> cam,
                       std::shared_ptr<Scene> scene, int spp, int max_depth)
    : camera(std::move(cam)), scene(std::move(scene)), spp(spp),
      max_depth(max_depth) {}

void print_ray(Ray &ray) {
	std::cout << "ray's origin is (" << ray.origin.x() << ", " << ray.origin.y()
	          << ", " << ray.origin.z() << ")" << std::endl;
	std::cout << "ray's direction is (" << ray.direction.x() << ", "
	          << ray.direction.y() << ", " << ray.direction.z() << ")"
	          << std::endl
	          << std::endl;
}

const Vec3f DUMMY_POS = Vec3f(0.f, 0.f, 0.f);

void Integrator::render() const {
	Vec2i resolution = camera->getImage()->getResolution();
	int cnt = 0;
	Sampler sampler;
	int sppx = (int) std::sqrt(spp);
	int sppy = sppx;
	auto start = std::chrono::steady_clock::now();
#ifndef MY_DEBUG
#pragma omp parallel for default(none), schedule(dynamic), \
        shared(cnt, resolution, start, std::cout), private(sampler)
#endif
	for (int dx = 0; dx < resolution.x(); dx++) {
#pragma omp critical
		{
			cnt++;
			UpdateProgress((float) cnt * 1.0f / (float) resolution.x());
			auto end = std::chrono::steady_clock::now();
			auto time =
			        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
			std::cout << ", Elapsed " << time << "s.";
			std::cout.flush();
		}

		sampler.setSeed(omp_get_thread_num());
		for (int dy = 0; dy < resolution.y(); dy++) {
			Vec3f L(0, 0, 0);
			// to compute radiance.
			for (int _spp = 0; _spp < spp; _spp++) {
				Vec2f cast_at = sampler.get2D() + Vec2f(dx, dy);
				Ray ray = camera->generateRay(cast_at.x(), cast_at.y());
				L += radiance(ray, sampler);
			}
			camera->getImage()->setPixel(dx, dy, L / spp);
		}
	}
}

Vec3f Integrator::radiance(Ray &ray, Sampler &sampler) const {
	Vec3f L(0, 0, 0);
	Vec3f beta(1, 1, 1);
	bool isDelta = false;
	Ray cur_ray = ray;
	for (int i = 0; i < max_depth; ++i) {
		/// Compute radiance (direct + indirect)
		Interaction interaction;
		scene->intersect(cur_ray, interaction);
		interaction.wo = -cur_ray.direction;
		if (interaction.type == Interaction::Type::NONE)
			break;
		else if (interaction.type == Interaction::Type::LIGHT) {
			if (i == 0) {
				L += scene->getLight()
				             ->emission(DUMMY_POS, interaction.wo)
				             .cwiseProduct(beta);
			}
			break;
		}

#ifdef MY_DEBUG
		std::cout << i << std::endl;
		DEBUG_VEC(interaction.normal);
		DEBUG_VEC(interaction.pos);
#endif
		// direct light
		L += directLighting(interaction, sampler).cwiseProduct(beta);

		// indirect light
		float pdf = interaction.material->sample(interaction, sampler);
		float cosine_wi = interaction.normal.dot(interaction.wi);
		Vec3f cur_brdf = interaction.material->evaluate(interaction);
		beta = beta.cwiseProduct(cur_brdf).eval() * cosine_wi / pdf;

		cur_ray = Ray(interaction.pos, interaction.wi);
	}
#ifdef MY_DEBUG
	DEBUG_VEC(L);
#endif
	return L;
}

Vec3f Integrator::directLighting(Interaction &interaction,
                                 Sampler &sampler) const {
	Vec3f L(0, 0, 0);
	// Compute direct lighting.
	float pdf = scene->getLight()->pdf(interaction, DUMMY_POS);
	Vec3f light_pos = scene->getLight()->sample(interaction, sampler);
	interaction.wi = (light_pos - interaction.pos).normalized();
	Ray shadow_ray(interaction.pos, interaction.wi);
	float dist = (light_pos - interaction.pos).norm();
	if (!scene->isShadowed(shadow_ray, dist)) {
		float dis_square = dist * dist;
		Vec3f cur_brdf = interaction.material->evaluate(interaction);
		float cosine_wi = interaction.normal.dot(interaction.wi);
		float cosine_theta = Vec3f(0, -1.0, 0).dot(-interaction.wi);
		L = scene->getLight()
		            ->emission(DUMMY_POS, -interaction.wi)
		            .cwiseProduct(cur_brdf) *
		    cosine_wi * cosine_theta / (pdf * dis_square);
	}

	return L;
}