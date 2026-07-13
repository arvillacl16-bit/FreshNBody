// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gravity.hpp"
#include "particle.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>

namespace fnb::accel {
  void basic(ParticleStore& particles, double epsilon2) {
    const size_t N = particles.N();
    const size_t test_thres = particles.test_thres().value_or(N);

    auto& positions = particles.positions;
    auto& accelerations = particles.accelerations;
    auto& mus = particles.mus;

    for (size_t i = 0; i < test_thres; ++i) {
      for (size_t j = i + 1; j < test_thres; ++j) {
        Vec3 dp = positions[i] - positions[j];
        double inv_r3 = 1. / (dp.mag2() + epsilon2);
        inv_r3 *= std::sqrt(inv_r3);
        accelerations[j] += dp * (mus[i] * inv_r3);
        accelerations[i] -= dp * (mus[j] * inv_r3);
      }

      for (size_t j = test_thres; j < N; ++j) {
        Vec3 dp = positions[i] - positions[j];
        double inv_r3 = 1. / (dp.mag2() + epsilon2);
        inv_r3 *= std::sqrt(inv_r3);
        accelerations[j] += dp * (mus[i] * inv_r3);
      }
    }
  }

  void compensated(ParticleStore& particles, double epsilon2) {
    const size_t N = particles.N();
    const size_t test_thres = particles.test_thres().value_or(N);

    auto& positions = particles.positions;
    auto& accelerations = particles.accelerations;
    auto& mus = particles.mus;

    for (size_t i = 0; i < N; ++i) {
      Vec3 sum;
      Vec3 c;
      for (size_t j = 0; j < test_thres; ++j) {
        if (i == j) continue;

        Vec3 dp = positions[i] - positions[j];
        double inv_r3 = 1. / (dp.mag2() + epsilon2);
        inv_r3 *= std::sqrt(inv_r3);
        Vec3 opp_force = dp * (mus[j] * inv_r3);

        Vec3 y = opp_force - c;
        Vec3 t = sum + y;
        c = (t - sum) - y;
        sum = t;
      }

      accelerations[i] = -sum;
    }
  }

  void jacobi(ParticleStore& particles, double epsilon2) {
    //
  }
} // namespace fnb::accel
