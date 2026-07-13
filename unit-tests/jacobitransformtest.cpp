#include "particle.hpp"
#include "transform.hpp"
#include <gtest/gtest.h>
#include <random>
#include <ranges>

namespace rv = std::views;

TEST(JacobiTransform, ReversibilityTest) {
  using namespace fnb;

  std::mt19937 gen(std::random_device{}());
  std::uniform_real_distribution<> dist(-20, 20);

  auto gen_random_particle = [&](unsigned long long id) {
    return IndParticle{
        .pos = {dist(gen) * 1.496e8, dist(gen) * 1.496e8, dist(gen) * 1.496e8},
        .vel = {dist(gen), dist(gen), dist(gen)},
        .mu = 398600,
        .id = id,
    };
  };

  ParticleStore p;
  ParticleStore useless;
  ParticleStore orig;

  p.add_particle({.mu = 1.327e11, .id = 0});
  useless.add_particle(IndParticle{.mu = 1.327e11});

  for (auto k : rv::iota(1, 7)) {
    p.add_particle(gen_random_particle(k));
    useless.add_particle(IndParticle{.mu = 398600});
  }

  orig = p;

  transform::inertial_to_jacobi_pos(p, useless);
  transform::inertial_to_jacobi_vel(p, useless);
  transform::jacobi_to_inertial_pos(p, useless);
  transform::jacobi_to_inertial_vel(p, useless);

  Vec3 total_error;
  for (auto i : rv::iota(0, 7)) {
    total_error += p[i].pos() - orig[i].pos();
    total_error += p[i].vel() - orig[i].vel();
  }

  ASSERT_LE(total_error.mag2(), 2.);
}
