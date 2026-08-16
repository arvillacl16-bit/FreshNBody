// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "simulation.hpp"
#include "particle.hpp"
#include <iostream>
#include <memory>
#include <unordered_map>

namespace fnb {
  struct Simulation::Impl {
    std::unordered_map<uint64_t, size_t> id_map;
  };

  Simulation::Simulation(Simulation&& other) noexcept = default;
  Simulation& Simulation::operator=(Simulation&&) noexcept = default;

  Simulation::Simulation(const Simulation& other)
      : pimpl_(std::make_unique<Impl>()), t_(other.t_), particles_(other.particles_), config(other.config), dt(other.dt) {
    pimpl_->id_map = other.pimpl_->id_map;
  }

  Simulation& Simulation::operator=(const Simulation& other) {
    if (this == &other) return *this;

    pimpl_ = std::make_unique<Impl>();
    t_ = other.t_;
    particles_ = other.particles_;
    dt = other.dt;
    config = other.config;

    pimpl_->id_map = other.pimpl_->id_map;
    return *this;
  }

  using Status = Simulation::Status;
  static std::string_view status_to_str(Status status) {
    switch (status) {
    case Status::TERM: return "TERM";
    default: return "";
    }
  }

  void Simulation::step(double dt) {
    switch (config.type()) {
    case IntegratorType::LEAPFROG: config.get<Leapfrog>().step(particles_, dt); break;
    case IntegratorType::WHFAST: config.get<WHFast>().step(particles_, dt); break;
    case IntegratorType::IAS15: config.get<IAS15>().step(particles_); break;
    case IntegratorType::MERCURIUS: config.get<Mercurius>().step(particles_, dt); break;
    default: int _;
    }
  }

  void Simulation::integrate(double t, bool exact) {
    while (t_ < t) {
      double step_size = exact ? (dt < t - t_ ? dt : t - t_) : dt;
      step(step_size);
      t_ += step_size;

      Status status = heartbeat(*this);
      if (status != Status::OK) std::cerr << "ERROR: Integration terminated with status " << status_to_str(status) << '\n';
    }
  }

  std::pair<Vec3, Vec3> Simulation::com() {
    const size_t N = particles_.test_thres().value_or(particles_.N());
    double total_mass = 0;
#pragma omp parallel for
    for (size_t i = 0; i < N; ++i) total_mass += particles_.mus[i];

    Vec3 pos{};
    Vec3 vel{};
#pragma omp parallel for reduction(+ : pos, vel)
    for (size_t i = 0; i < particles_.test_thres().value_or(particles_.N()); ++i) {
      pos += particles_.positions[i] * particles_.mus[i] / total_mass;
      vel += particles_.velocities[i] * particles_.mus[i] / total_mass;
    }
    return {pos, vel};
  }

  void Simulation::move_to_com() {
    auto res = com();
#pragma omp parallel for
    for (size_t i = 0; i < particles_.N(); ++i) {
      particles_.positions[i] -= res.first;
      particles_.velocities[i] -= res.second;
    }
  }
} // namespace fnb
