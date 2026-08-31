#include <GTest/gtest.h>
#include "particle.hpp"
#include "integrator.hpp"
#include "gravity.hpp"

TEST(IntegratorTests, Ias15) {  
  using namespace fnb;
  ParticleStore sim;

  Vec3 correct_pos{1, 0, 0};
  sim.add_particles(
    {{.pos = {0, 0, 0}, .vel = {0, 0, 0}, .mu = 4 * M_PI * M_PI, .id = 0},
      {.pos = correct_pos, .vel = {0, 100, 0}, .mu = 0.001, .id = 1},
      {.pos = {1.5, 0, 0}, .vel = {0, 1, 0}, .mu = 0.0001, .id = 2}});

  IAS15 integ;
  integ.epsilon = 1e-6;
  for (size_t i = 0; i < 1000; ++i) {
    integ.step(sim);
  }

  ASSERT_LE(sim[1].pos().dist2(correct_pos), 1e14);

}
