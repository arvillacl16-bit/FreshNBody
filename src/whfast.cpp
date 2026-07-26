// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gravity.hpp"
#include "integrator.hpp"
#include "particle.hpp"
#include "transform.hpp"
#include <unistd.h>

namespace fnb {
  namespace correctors {
    constexpr double a_1 = 0.41833001326703777398908601289259374469640768464934;
    constexpr double a_2 = 0.83666002653407554797817202578518748939281536929867;
    constexpr double a_3 = 1.2549900398011133219672580386777812340892230539480;
    constexpr double a_4 = 1.6733200530681510959563440515703749787856307385973;
    constexpr double a_5 = 2.0916500663351888699454300644629687234820384232467;
    constexpr double a_6 = 2.5099800796022266439345160773555624681784461078960;
    constexpr double a_7 = 2.9283100928692644179236020902481562128748537925454;
    constexpr double a_8 = 3.3466401061363021919126881031407499575712614771947;
    constexpr double b_31 = -0.024900596027799867499350357910273437184309981229127;
    constexpr double b_51 = -0.0083001986759332891664501193034244790614366604097090;
    constexpr double b_52 = 0.041500993379666445832250596517122395307183302048545;
    constexpr double b_71 = 0.0024926811426922105779030593952776964450539008582219;
    constexpr double b_72 = -0.018270923246702131478062356884535264841652263842597;
    constexpr double b_73 = 0.053964399093127498721765893493510877532452806339655;
    constexpr double b_111 = 0.00020361579647854651301632818774633716473696537436847;
    constexpr double b_112 = -0.0023487215292295354188307328851055489876255097419754;
    constexpr double b_113 = 0.012309078592019946317544564763237909911330686448336;
    constexpr double b_114 = -0.038121613681288650508647613260247372125243616270670;
    constexpr double b_178 = 0.093056103771425958591541059067553547100903397724386;
    constexpr double b_177 = -0.065192863576377893658290760803725762027864651086787;
    constexpr double b_176 = 0.032422198864713580293681523029577130832258806467604;
    constexpr double b_175 = -0.012071760822342291062449751726959664253913904872527;
    constexpr double b_174 = 0.0033132577069380655655490196833451994080066801611459;
    constexpr double b_173 = -0.00063599983075817658983166881625078545864140848560259;
    constexpr double b_172 = 0.000076436355227935738363241846979413475106795392377415;
    constexpr double b_171 = -0.0000043347415473373580190650223498124944896789841432241;
    constexpr double corrector_2 = 0.03486083443891981449909050107438281205803;
  } // namespace correctors

  static constexpr double invfactorial[35] = {
      1.,
      1.,
      1. / 2.,
      1. / 6.,
      1. / 24.,
      1. / 120.,
      1. / 720.,
      1. / 5040.,
      1. / 40320.,
      1. / 362880.,
      1. / 3628800.,
      1. / 39916800.,
      1. / 479001600.,
      1. / 6227020800.,
      1. / 87178291200.,
      1. / 1307674368000.,
      1. / 20922789888000.,
      1. / 355687428096000.,
      1. / 6402373705728000.,
      1. / 121645100408832000.,
      1. / 2432902008176640000.,
      1. / 51090942171709440000.,
      1. / 1124000727777607680000.,
      1. / 25852016738884976640000.,
      1. / 620448401733239439360000.,
      1. / 15511210043330985984000000.,
      1. / 403291461126605635584000000.,
      1. / 10888869450418352160768000000.,
      1. / 304888344611713860501504000000.,
      1. / 8841761993739701954543616000000.,
      1. / 265252859812191058636308480000000.,
      1. / 8222838654177922817725562880000000.,
      1. / 263130836933693530167218012160000000.,
      1. / 8683317618811886495518194401280000000.,
      1. / 295232799039604140847618609643520000000.};

  namespace {
    template <typename T> T max(T a, T b) { return a < b ? b : a; }
    template <typename T> T min(T a, T b) { return a < b ? a : b; }

    constexpr double fastabs(double x) { return x > 0 ? x : -x; }
    void stumpff_cs(double* restrict cs, double z) {
      unsigned int n = 0;
      while (fastabs(z) > 0.1) {
        z = z / 4.;
        ++n;
      }

      constexpr int nmax = 15;
      double c_odd = invfactorial[nmax];
      double c_even = invfactorial[nmax - 1];

      for (int np = nmax - 2; np >= 5; np -= 2) {
        c_odd = invfactorial[np] - z * c_odd;
        c_even = invfactorial[np - 1] - z * c_even;
      }

      cs[5] = c_odd;
      cs[4] = c_even;
      cs[3] = invfactorial[3] - z * cs[5];
      cs[2] = invfactorial[2] - z * cs[4];
      cs[1] = invfactorial[1] - z * cs[3];
      for (; n > 0; n--) {
        z *= 4.;
        cs[5] = (cs[5] + cs[4] + cs[3] * cs[2]) * 0.0625;
        cs[4] = (1. + cs[1]) * cs[3] * 0.125;
        cs[3] = 1. / 6. - z * cs[5];
        cs[2] = 0.5 - z * cs[4];
        cs[1] = 1. - z * cs[3];
      }
      cs[0] = invfactorial[0] - z * cs[2];
    }

    void stumpff_cs3(double* restrict cs, double z) {
      unsigned int n = 0;
      while (fastabs(z) > 0.1) {
        z /= 4.;
        ++n;
      }

      constexpr int nmax = 13;
      double c_odd = invfactorial[nmax];
      double c_even = invfactorial[nmax - 1];
      for (int np = nmax - 2; np >= 3; np -= 2) {
        c_odd = invfactorial[np] - z * c_odd;
        c_even = invfactorial[np - 1] - z * c_even;
      }

      cs[3] = c_odd;
      cs[2] = c_even;
      cs[1] = invfactorial[1] - z * c_odd;
      cs[0] = invfactorial[0] - z * c_even;
      for (; n > 0; --n) {
        cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
        cs[2] = cs[1] * cs[1] * 0.5;
        cs[1] = cs[0] * cs[1];
        cs[0] = 2. * cs[0] * cs[0] - 1.;
      }
    }

    void stiefel_Gs(double* restrict Gs, double beta, double X) {
      double X2 = X * X;
      stumpff_cs(Gs, beta * X2);
      Gs[1] *= X;
      Gs[2] *= X2;
      double pow = X2 * X;
      Gs[3] *= pow;
      pow *= X;
      Gs[4] *= pow;
      pow *= X;
      Gs[5] *= pow;
      return;
    }

    void stiefel_Gs3(double* restrict Gs, double beta, double X) {
      double X2 = X * X;
      stumpff_cs3(Gs, beta * X2);
      Gs[1] *= X;
      Gs[2] *= X2;
      Gs[3] *= X * X2;
      return;
    }

    void update_accel(GravityMethod method, double epsilon2, ParticleStore& particles) {
      switch (method) {
      case GravityMethod::BASIC: accel::basic(particles, epsilon2); break;
      case GravityMethod::COMPENSATED: accel::compensated(particles, epsilon2); break;
      case GravityMethod::JACOBI: accel::jacobi(particles, epsilon2); break;
      default: break;
      }
    }

    void init_dummy_store(WHFast& obj, ParticleStore& particles) {
      while (obj.p_j.N() < particles.N()) {
        obj.p_j.add_particle(IndParticle{});
      }
    }
  } // namespace

  void WHFast::step(ParticleStore& particles, double dt) {
    init_dummy_store(*this, particles);

    // Helper lambda performing a single WH substep with a scaled time interval
    auto wh_substep = [&](double sub_dt) {
      if (sub_dt == 0.0) return;

      update_accel(gravity, epsilon * epsilon, particles);
      double star_mu = particles.mus[0];
      Vec3 star_pos = particles.positions[0];

#pragma omp parallel for
      for (size_t i = 1; i < particles.N(); ++i) {
        Vec3 dr = star_pos - particles.positions[i];
        particles.velocities[i] += sub_dt / 2. * (particles.accelerations[i] - star_mu * dr / dr.mag());
      }

      transform::inertial_to_jacobi_pos(particles, p_j);
      transform::inertial_to_jacobi_vel(particles, p_j);

      double interior_mu = star_mu;
      for (size_t i = 1; i < particles.N(); ++i) {
        double r = p_j.positions[i].mag();
        double vr = p_j.positions[i].dot(p_j.velocities[i]) / r;
        double alpha = 2 / r - p_j.velocities[i].mag2() / interior_mu;

        double chi = sqrt(interior_mu) * sub_dt * alpha;
        double sqrt_mu = sqrt(interior_mu);

        for (int j = 0; j < 10; ++j) {
          double z = alpha * chi * chi;
          double cs[6];
          stumpff_cs(cs, z);

          const double C = cs[0];
          const double S = cs[1];

          double f = (r * vr / sqrt_mu) * chi * chi * C + (1.0 - alpha * r) * chi * chi * chi * S + r * chi - sqrt_mu * sub_dt;
          double df = (r * vr / sqrt_mu) * chi * (1.0 - alpha * chi * chi * S) + (1.0 - alpha * r) * chi * chi * C + r;

          double dchi = f / df;
          chi -= dchi;

          if (fastabs(dchi) < 1e-12) break;
        }

        double cs[6];
        stumpff_cs(cs, alpha * chi * chi);
        const double C = cs[0];
        const double S = cs[1];

        double f = 1 - chi * chi / r * C;
        double g = sub_dt - chi * chi * chi / sqrt_mu * S;

        Vec3 old_p = p_j.positions[i];
        Vec3 old_v = p_j.velocities[i];
        p_j.positions[i] = f * old_p + g * old_v;

        double r_new = p_j.positions[i].mag();
        double f_dot = (sqrt_mu / (r * r_new)) * (alpha * chi * chi * chi * S - chi);
        double g_dot = 1.0 - (chi * chi / r_new) * C;

        p_j.velocities[i] = f_dot * old_p + g_dot * old_v;

        interior_mu += particles.mus[i];
      }

      transform::jacobi_to_inertial_vel(p_j, particles);
      transform::jacobi_to_inertial_pos(p_j, particles);

      star_pos = particles.positions[0];
#pragma omp parallel for
      for (size_t i = 1; i < particles.N(); ++i) {
        Vec3 dr = star_pos - particles.positions[i];
        particles.velocities[i] += sub_dt / 2. * (particles.accelerations[i] - star_mu * dr / dr.mag());
      }
    };

    auto apply_corrector_kick = [&](double b_coeff) {
      if (b_coeff == 0.0) return;
      double star_mu = particles.mus[0];
      Vec3 star_pos = particles.positions[0];
#pragma omp parallel for
      for (size_t i = 1; i < particles.N(); ++i) {
        Vec3 dr = star_pos - particles.positions[i];
        particles.velocities[i] += b_coeff * dt * (particles.accelerations[i] - star_mu * dr / dr.mag());
      }
    };

    using namespace correctors;

    wh_substep(a_1 * dt);
    apply_corrector_kick(b_31);
    wh_substep((a_2 - a_1) * dt);
    wh_substep((a_3 - a_2) * dt);
    apply_corrector_kick(b_51);
    wh_substep((a_4 - a_3) * dt);
    apply_corrector_kick(b_52);
    wh_substep((a_5 - a_4) * dt);
    apply_corrector_kick(b_71);
    wh_substep((a_6 - a_5) * dt);
    apply_corrector_kick(b_72);
    wh_substep((a_7 - a_6) * dt);
    apply_corrector_kick(b_73);
    wh_substep((a_8 - a_7) * dt);
    wh_substep((1.0 - a_8) * dt);
  }
} // namespace fnb
