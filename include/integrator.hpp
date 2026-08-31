// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

#include "particle.hpp"
#include <cstdint>

namespace fnb {
  enum class GravityMethod : uint8_t { NONE, BASIC, JACOBI, COMPENSATED };

  class Leapfrog {
  public:
    GravityMethod gravity = GravityMethod::BASIC;
    double epsilon = 1e-12;

    void step(ParticleStore& particles, double dt);
  };

  class WHFast {
  private:
    ParticleStore p_j;
  public:
    GravityMethod gravity = GravityMethod::BASIC;
    double epsilon = 1e-12;


    void step(ParticleStore& particles, double dt);
  };

  class IAS15 {
  private:
    double dt_ = 10000;
    std::vector<Vec3> y0_, v0_, a0_;
    std::vector<Vec3> a_sub_, g_;
    std::vector<Vec3> current_accs_;
  public:
    unsigned int max_iterations = 12;
    double epsilon = 1e-12;
    double epsilon_dtadjust = 1e-9;
    void step(ParticleStore& particles);
  };

  class Mercurius {
  public:
    void step(ParticleStore& particles, double dt);
  };
} // namespace fnb

#endif
