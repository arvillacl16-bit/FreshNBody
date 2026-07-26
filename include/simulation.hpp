// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include "integrator.hpp"
#include "particle.hpp"
#include <cstdint>
#include <functional>
#include <optional>

namespace fnb {
  enum class IntegratorType : uint8_t { NONE = 0, LEAPFROG, WHFAST, IAS15, MERCURIUS };

  struct NoIntegration {};

  class BadConfigAccess : public std::exception {
  public:
    const char* what() const noexcept override { return "Wrong type used"; }
  };

  class Simulation {
  private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;

    double t_;

    ParticleStore particles_;

  public:
    class IntegratorConfig {
    private:
      IntegratorType type_;
      union {
        Leapfrog leapfrog_;
        WHFast whfast_;
        IAS15 ias15_;
        Mercurius mercurius_;
        NoIntegration nointeg_;
      };

    public:
      IntegratorConfig() : type_(IntegratorType::NONE), nointeg_{} {}
      explicit IntegratorConfig(IntegratorType type) : type_(type) { reset(); }

      IntegratorConfig(const Leapfrog& leapfrog) : type_(IntegratorType::LEAPFROG), leapfrog_(leapfrog) {}
      IntegratorConfig(const WHFast& whfast) : type_(IntegratorType::WHFAST), whfast_(whfast) {}
      IntegratorConfig(const IAS15& ias15) : type_(IntegratorType::IAS15), ias15_(ias15) {}
      IntegratorConfig(const Mercurius& mercurius) : type_(IntegratorType::MERCURIUS), mercurius_(mercurius) {}

      IntegratorConfig(IntegratorConfig&& other) noexcept {
        switch (other.type_) {
        case IntegratorType::NONE: nointeg_ = other.nointeg_; break;
        case IntegratorType::LEAPFROG: leapfrog_ = other.leapfrog_; break;
        case IntegratorType::WHFAST: whfast_ = other.whfast_; break;
        case IntegratorType::IAS15: ias15_ = other.ias15_; break;
        case IntegratorType::MERCURIUS: mercurius_ = other.mercurius_; break;
        }
        type_ = other.type_;
        other.reset();
      }

      IntegratorConfig& operator=(IntegratorConfig&& other) noexcept {
        if (this == &other) return *this;
        switch (other.type_) {
        case IntegratorType::NONE: nointeg_ = other.nointeg_; break;
        case IntegratorType::LEAPFROG: leapfrog_ = other.leapfrog_; break;
        case IntegratorType::WHFAST: whfast_ = other.whfast_; break;
        case IntegratorType::IAS15: ias15_ = other.ias15_; break;
        case IntegratorType::MERCURIUS: mercurius_ = other.mercurius_; break;
        }
        type_ = other.type_;
        other.reset();
        return *this;
      }

      IntegratorConfig(const IntegratorConfig& other) {
        switch (other.type_) {
        case IntegratorType::NONE: nointeg_ = other.nointeg_; break;
        case IntegratorType::LEAPFROG: leapfrog_ = other.leapfrog_; break;
        case IntegratorType::WHFAST: whfast_ = other.whfast_; break;
        case IntegratorType::IAS15: ias15_ = other.ias15_; break;
        case IntegratorType::MERCURIUS: mercurius_ = other.mercurius_; break;
        }
        type_ = other.type_;
      }

      IntegratorConfig& operator=(const IntegratorConfig& other) noexcept {
        if (this == &other) return *this;
        switch (other.type_) {
        case IntegratorType::NONE: nointeg_ = other.nointeg_; break;
        case IntegratorType::LEAPFROG: leapfrog_ = other.leapfrog_; break;
        case IntegratorType::WHFAST: whfast_ = other.whfast_; break;
        case IntegratorType::IAS15: ias15_ = other.ias15_; break;
        case IntegratorType::MERCURIUS: mercurius_ = other.mercurius_; break;
        }
        type_ = other.type_;
        return *this;
      }

      ~IntegratorConfig() {
        switch (type_) {
        case IntegratorType::WHFAST: whfast_.~WHFast();
        default: break;
        }
      }

      void set_type(IntegratorType type, bool do_reset = true) {
        type_ = type;
        if (do_reset) reset();
      }

      IntegratorType type() { return type_; }

      void reset() {
        switch (type_) {
        case IntegratorType::LEAPFROG: leapfrog_ = Leapfrog{}; break;
        case IntegratorType::WHFAST: whfast_ = WHFast{}; break;
        case IntegratorType::IAS15: ias15_ = IAS15{}; break;
        case IntegratorType::MERCURIUS: mercurius_ = Mercurius{}; break;
        default: nointeg_ = NoIntegration{};
        }
      }

      template <typename Integrator> Integrator& get() {
        if constexpr (std::is_same_v<Integrator, Leapfrog>) {
          if (type_ != IntegratorType::LEAPFROG) throw BadConfigAccess{};
          return leapfrog_;
        } else if constexpr (std::is_same_v<Integrator, WHFast>) {
          if (type_ != IntegratorType::WHFAST) throw BadConfigAccess{};
          return whfast_;
        } else if constexpr (std::is_same_v<Integrator, IAS15>) {
          if (type_ != IntegratorType::IAS15) throw BadConfigAccess{};
          return ias15_;
        } else if constexpr (std::is_same_v<Integrator, Mercurius>) {
          if (type_ != IntegratorType::MERCURIUS) throw BadConfigAccess{};
          return mercurius_;
        } else
          static_assert(sizeof(Integrator) == 0, "Invalid type");
      }

      template <typename Integrator> const Integrator& get() const {
        if constexpr (std::is_same_v<Integrator, Leapfrog>) {
          if (type_ != IntegratorType::LEAPFROG) throw BadConfigAccess{};
          return leapfrog_;
        } else if constexpr (std::is_same_v<Integrator, WHFast>) {
          if (type_ != IntegratorType::WHFAST) throw BadConfigAccess{};
          return whfast_;
        } else if constexpr (std::is_same_v<Integrator, IAS15>) {
          if (type_ != IntegratorType::IAS15) throw BadConfigAccess{};
          return ias15_;
        } else if constexpr (std::is_same_v<Integrator, Mercurius>) {
          if (type_ != IntegratorType::MERCURIUS) throw BadConfigAccess{};
          return mercurius_;
        } else
          static_assert(sizeof(Integrator) == 0, "Invalid type");
      }
    } config;

    enum class Status : uint8_t { OK = 0, TERM };

    double dt;

    std::function<Status(Simulation&)> heartbeat;

    Simulation() = delete;
    Simulation(const Simulation&);
    Simulation(Simulation&&) noexcept;
    Simulation& operator=(const Simulation&);
    Simulation& operator=(Simulation&&) noexcept;
    ~Simulation();

    void step(double dt);
    void steps(size_t n, std::optional<double> dt_ = std::optional<double>{}) {
      double step_size = dt_.value_or(dt);
      for (size_t i = 0; i < n; ++i) step(step_size);
    }

    void integrate(double t, bool exact = false);

    void add_particle(const IndParticle& p) { particles_.add_particle(p); }
    void add_particles(const std::vector<IndParticle>& ps) { particles_.add_particles(ps); }
    void remove_particle(size_t idx) { particles_.remove_particle(idx); }

    std::pair<Vec3, Vec3> com();
    void move_to_com();
  };
} // namespace fnb

#endif
