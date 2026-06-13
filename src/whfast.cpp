// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrator.hpp"
#include "particle.hpp"

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
  }
  
  static constexpr double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

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
      for(;n>0;n--) {
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
      for (;n > 0; --n) {
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
  }

  
  void WHFast::step(ParticleStore& particles, double dt) {
    //
  }
} // namespace fnb
