// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "transform.hpp"
#include "particle.hpp"
#include <cstddef>

namespace fnb::transform {
  void inertial_to_jacobi_pos(const ParticleStore& restrict particles, ParticleStore& restrict p_j) {
    const size_t N = particles.N();
    const size_t N_thres = particles.test_thres().value_or(N);
    if (N < 2) return;

    Vec3 current_com_pos = particles.mus[0] * particles.positions[0];
    double current_tot_mass = particles.mus[0];
    for (size_t i = 1; i < N_thres; ++i) {
      p_j.positions[i] = particles.positions[i] - current_com_pos / current_tot_mass;

      current_com_pos += particles.mus[i] * particles.positions[i];
      current_tot_mass += particles.mus[i];
    }

    if (N_thres < N) {
      current_com_pos += particles.mus[N_thres] * particles.positions[N_thres];
      current_tot_mass += particles.mus[N_thres];
    }

    Vec3 act_com_pos = current_com_pos / current_tot_mass;
    for (size_t i = N_thres; i < N; ++i) {
      p_j.positions[i] = particles.positions[i] - act_com_pos;
    }

    p_j.positions[0] = act_com_pos;
  }

  void inertial_to_jacobi_vel(const ParticleStore& restrict particles, ParticleStore& restrict p_j) {
    const size_t N = particles.N();
    const size_t N_thres = particles.test_thres().value_or(N);
    if (N < 2) return;

    Vec3 current_com_pos = particles.mus[0] * particles.positions[0];
    double current_tot_mass = particles.mus[0];
    for (size_t i = 1; i < N_thres; ++i) {
      p_j.velocities[i] = particles.velocities[i] - current_com_pos / current_tot_mass;

      current_com_pos += particles.mus[i] * particles.velocities[i];
      current_tot_mass += particles.mus[i];
    }

    if (N_thres < N) {
      current_com_pos += particles.mus[N_thres] * particles.velocities[N_thres];
      current_tot_mass += particles.mus[N_thres];
    }

    Vec3 act_com_pos = current_com_pos / current_tot_mass;
    for (size_t i = N_thres; i < N; ++i) {
      p_j.velocities[i] = particles.velocities[i] - act_com_pos;
    }

    p_j.velocities[0] = act_com_pos;
  }

  void inertial_to_jacobi_acc(const ParticleStore& restrict particles, ParticleStore& restrict p_j) {
    const size_t N = particles.N();
    const size_t N_thres = particles.test_thres().value_or(N);
    if (N < 2) return;

    Vec3 current_com_pos = particles.mus[0] * particles.positions[0];
    double current_tot_mass = particles.mus[0];
    for (size_t i = 1; i < N_thres; ++i) {
      p_j.accelerations[i] = particles.accelerations[i] - current_com_pos / current_tot_mass;

      current_com_pos += particles.mus[i] * particles.accelerations[i];
      current_tot_mass += particles.mus[i];
    }

    if (N_thres < N) {
      current_com_pos += particles.mus[N_thres] * particles.accelerations[N_thres];
      current_tot_mass += particles.mus[N_thres];
    }

    Vec3 act_com_pos = current_com_pos / current_tot_mass;
    for (size_t i = N_thres; i < N; ++i) {
      p_j.accelerations[i] = particles.accelerations[i] - act_com_pos;
    }

    p_j.accelerations[0] = act_com_pos;
  }

  void jacobi_to_inertial_pos(ParticleStore& restrict particles, const ParticleStore& restrict p_j) {
    const size_t N = particles.N();
    if (N < 2) return;

    const size_t N_thres = particles.test_thres().value_or(N);
    const Vec3 act_com_pos = p_j.positions[0];

    double total_chain_mass = p_j.mus[0];
    for (size_t i = 1; i < N_thres; ++i) {
      total_chain_mass += p_j.mus[i];
    }
    if (N_thres < N) {
      total_chain_mass += p_j.mus[N_thres];
    }

    double current_tot_mass = total_chain_mass;
    Vec3 current_com_pos = act_com_pos;

    for (size_t i = N_thres - 1; i > 0; --i) {
      double sub_mass = current_tot_mass - p_j.mus[i];

      particles.positions[i] = current_com_pos + (sub_mass / current_tot_mass) * p_j.positions[i];

      current_com_pos = current_com_pos - (p_j.mus[i] / current_tot_mass) * p_j.positions[i];
      current_tot_mass = sub_mass;
    }

    particles.positions[0] = current_com_pos;

#pragma omp parallel for schedule(static) if (N - N_thres > 64)
    for (size_t i = N_thres; i < N; ++i) {
      particles.positions[i] = p_j.positions[i] + act_com_pos;
    }
  }

  void jacobi_to_inertial_vel(ParticleStore& restrict particles, const ParticleStore& restrict p_j) {
    const size_t N = particles.N();
    if (N < 2) return;

    const size_t N_thres = particles.test_thres().value_or(N);
    const Vec3 act_com_pos = p_j.velocities[0];

    double total_chain_mass = p_j.mus[0];
    for (size_t i = 1; i < N_thres; ++i) {
      total_chain_mass += p_j.mus[i];
    }
    if (N_thres < N) {
      total_chain_mass += p_j.mus[N_thres];
    }

    double current_tot_mass = total_chain_mass;
    Vec3 current_com_pos = act_com_pos;

    for (size_t i = N_thres - 1; i > 0; --i) {
      double sub_mass = current_tot_mass - p_j.mus[i];

      particles.velocities[i] = current_com_pos + (sub_mass / current_tot_mass) * p_j.velocities[i];

      current_com_pos = current_com_pos - (p_j.mus[i] / current_tot_mass) * p_j.velocities[i];
      current_tot_mass = sub_mass;
    }

    particles.velocities[0] = current_com_pos;

#pragma omp parallel for schedule(static) if (N - N_thres > 64)
    for (size_t i = N_thres; i < N; ++i) {
      particles.velocities[i] = p_j.velocities[i] + act_com_pos;
    }
  }

  void inertial_to_demhelio_pos(const ParticleStore& restrict particles, ParticleStore& restrict p_h) {
    double total_mu = 0;
    Vec3 bary;
#pragma omp parallel for reduction(+ : bary, total_mu)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.mus[i] * particles.positions[i];
      total_mu += particles.mus[i];
    }
    p_h.positions[0] = bary / total_mu;
    Vec3 star_pos = particles.positions[0];
#pragma omp parallel for default(none) shared(p_h, particles, star_pos)
    for (size_t i = 1; i < particles.N(); ++i) p_h.positions[i] = particles.positions[i] - star_pos;
  }

  void inertial_to_demhelio_vel(const ParticleStore& restrict particles, ParticleStore& restrict p_h) {
    double total_mu = 0;
    Vec3 bary;
#pragma omp parallel for reduction(+ : bary, total_mu)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.mus[i] * particles.velocities[i];
      total_mu += particles.mus[i];
    }
    p_h.velocities[0] = bary / total_mu;
    Vec3 star_pos = particles.velocities[0];
#pragma omp parallel for default(none) shared(p_h, particles, star_pos)
    for (size_t i = 1; i < particles.N(); ++i) p_h.velocities[i] = particles.velocities[i] - star_pos;
  }

  void inertial_to_demhelio_acc(const ParticleStore& restrict particles, ParticleStore& restrict p_h) {
    double total_mu = 0;
    Vec3 bary;
#pragma omp parallel for reduction(+ : bary, total_mu)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.mus[i] * particles.accelerations[i];
      total_mu += particles.mus[i];
    }
    p_h.accelerations[0] = bary / total_mu;
    Vec3 star_pos = particles.accelerations[0];
#pragma omp parallel for default(none) shared(p_h, particles, star_pos)
    for (size_t i = 1; i < particles.N(); ++i) p_h.accelerations[i] = particles.accelerations[i] - star_pos;
  }

  void demhelio_to_inertial_pos(ParticleStore& restrict particles, const ParticleStore& restrict p_h) {
    double total_mu = particles.mus[0];
    Vec3 everything_except_star;
#pragma omp parallel for reduction(+ : total_mu, everything_except_star)
    for (size_t i = 1; i < particles.test_thres().value_or(particles.N()); ++i) {
      everything_except_star += p_h.positions[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }

    Vec3 x0 = particles.positions[0] = p_h.positions[0] - everything_except_star / total_mu;
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) particles.positions[i] = p_h.positions[i] + x0;
  }

  void demhelio_to_inertial_vel(ParticleStore& restrict particles, const ParticleStore& restrict p_h) {
    double total_mu = particles.mus[0];
    Vec3 everything_except_star;
#pragma omp parallel for reduction(+ : total_mu, everything_except_star)
    for (size_t i = 1; i < particles.test_thres().value_or(particles.N()); ++i) {
      everything_except_star += p_h.velocities[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }

    Vec3 x0 = particles.velocities[0] = p_h.velocities[0] - everything_except_star / total_mu;
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) particles.velocities[i] = p_h.velocities[i] + x0;
  }

  void demhelio_to_inertial_acc(ParticleStore& restrict particles, const ParticleStore& restrict p_h) {
    double total_mu = particles.mus[0];
    Vec3 everything_except_star;
#pragma omp parallel for reduction(+ : total_mu, everything_except_star)
    for (size_t i = 1; i < particles.test_thres().value_or(particles.N()); ++i) {
      everything_except_star += p_h.accelerations[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }

    Vec3 x0 = particles.accelerations[0] = p_h.accelerations[0] - everything_except_star / total_mu;
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) particles.accelerations[i] = p_h.accelerations[i] + x0;
  }

  void inertial_to_whds_pos(const ParticleStore& restrict particles, ParticleStore& restrict p_h) {
    Vec3 star_pos = p_h.positions[0] = particles.positions[0];
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) p_h.positions[i] = particles.positions[i] - star_pos;
  }

  void inertial_to_whds_vel(const ParticleStore& restrict particles, ParticleStore& restrict p_h) {
    Vec3 total_momentum;
    double total_mu = 0;
#pragma omp parallel for reduction(+ : total_momentum, total_mu)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      total_momentum += particles.mus[i] * particles.velocities[i];
      total_mu += particles.mus[i];
    }

    p_h.velocities[0] = total_momentum / total_mu;
    Vec3 star_vel = particles.velocities[0];
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) p_h.velocities[i] = particles.velocities[i] - star_vel;
  }

  void whds_to_inertial_pos(ParticleStore& restrict particles, const ParticleStore& restrict p_h) {
    Vec3 star_pos = particles.positions[0] = p_h.positions[0];
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) particles.positions[i] = star_pos + p_h.positions[i];
  }

  void whds_to_inertial_vel(ParticleStore& restrict particles, const ParticleStore& restrict p_h) {
    Vec3 star_vel_correction;
    double total_mu = 0;
#pragma omp parallel for reduction(+ : star_vel_correction, total_mu)
    for (size_t i = 1; i < particles.N(); ++i) {
      star_vel_correction += p_h.mus[i] * p_h.velocities[i];
      total_mu += p_h.mus[i];
    }
    total_mu += particles.mus[0];

    Vec3 v0 = particles.velocities[0] = p_h.velocities[0] - star_vel_correction / total_mu;
#pragma omp parallel for
    for (size_t i = 1; i < particles.N(); ++i) particles.velocities[i] = v0 + p_h.velocities[i];
  }

  void inertial_to_barycentric_pos(const ParticleStore& restrict particles, ParticleStore& restrict p_b) {
    Vec3 bary;
    double total_mu = 0;
#pragma omp parallel for reduction(+ : bary, total_mu) if (particles.N() > 64)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.positions[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }
    bary /= total_mu;
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) p_b.positions[i] = particles.positions[i] - bary;
  }

  void inertial_to_barycentric_vel(const ParticleStore& restrict particles, ParticleStore& restrict p_b) {
    Vec3 bary;
    double total_mu = 0;
#pragma omp parallel for reduction(+ : bary, total_mu) if (particles.N() > 64)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.velocities[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }
    bary /= total_mu;
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) p_b.velocities[i] = particles.velocities[i] - bary;
  }

  void inertial_to_barycentric_acc(const ParticleStore& restrict particles, ParticleStore& restrict p_b) {
    Vec3 bary;
    double total_mu = 0;
#pragma omp parallel for reduction(+ : bary, total_mu) if (particles.N() > 64)
    for (size_t i = 0; i < particles.test_thres().value_or(particles.N()); ++i) {
      bary += particles.accelerations[i] * particles.mus[i];
      total_mu += particles.mus[i];
    }
    bary /= total_mu;
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) p_b.accelerations[i] = particles.accelerations[i] - bary;
  }

  void barycentric_to_inertial_pos(ParticleStore& restrict particles, const ParticleStore& restrict p_b) {
    Vec3 star = p_b.positions[0];
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) particles.positions[i] = p_b.positions[i] - star;
  }

  void barycentric_to_inertial_vel(ParticleStore& restrict particles, const ParticleStore& restrict p_b) {
    Vec3 star = p_b.velocities[0];
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) particles.velocities[i] = p_b.velocities[i] - star;
  }

  void barycentric_to_inertial_acc(ParticleStore& restrict particles, const ParticleStore& restrict p_b) {
    Vec3 star = p_b.accelerations[0];
#pragma omp parallel for if (particles.N() > 64)
    for (size_t i = 0; i < particles.N(); ++i) particles.accelerations[i] = p_b.accelerations[i] - star;
  }
} // namespace fnb::transform
