#include "integrator.hpp"
#include "particle.hpp"
#include <cmath>
#include <gtest/gtest.h>

constexpr double square(double r) { return r * r; }

TEST(Logic, LeapfrogBasic) {
  using namespace fnb;
  ParticleStore sim;

  Vec3 correct_pos{1, 0, 0};
  sim.add_particles(
      {{.pos = {0, 0, 0}, .vel = {0, 0, 0}, .mu = 4 * M_PI * M_PI, .id = 0},
       {.pos = correct_pos, .vel = {0, 1, 0}, .mu = 0.001, .id = 1},
       {.pos = {1.5, 0, 0}, .vel = {0, 1, 0}, .mu = 0.0001, .id = 2}});

  WHFast integ;
  integ.gravity = GravityMethod::BASIC;
  integ.epsilon = 1e-6;
  constexpr double step_size = 1. / 40;
  for (int k = 0; k < 40; ++k) {
    integ.step(sim, step_size);
  }

  ASSERT_LE(sim[1].pos().dist2(correct_pos), square(1e7));
}
