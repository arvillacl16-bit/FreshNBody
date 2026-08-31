#include "integrator.hpp"
#include "particle.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

namespace fnb {
  namespace {
    static constexpr double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626 };
    static constexpr double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147};
    static constexpr double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};
    static constexpr double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588};
    static constexpr double w[8] = {0.03125, 0.185358154802979278540728972807180754479812609, 0.304130620646785128975743291458180383736715043, 0.376517545389118556572129261157225608762708603, 0.391572167452493593082499533303669362149363727, 0.347014795634501068709955597003528601733139176, 0.249647901329864963257869294715235590174262844, 0.114508814744257199342353731044292225247093225};

    void compute_accelerations(ParticleStore& particles, std::vector<Vec3>& accs) {
      size_t N = particles.N();
      std::fill(accs.begin(), accs.end(), Vec3{0, 0, 0});
      for (size_t i = 0; i < N; ++i) {
        Vec3 pi = particles.positions[i];
        for (size_t j = i + 1; j < N; ++j) {
          Vec3 pj = particles.positions[j];
          Vec3 rij = pj - pi;
          double dist2 = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
          double dist = std::sqrt(dist2);
          if (dist < 1e-9) continue;
          double inv_dist3 = 1.0 / (dist2 * dist);
          accs[i] += rij * (particles.mus[j] * inv_dist3);
          accs[j] += rij * (-particles.mus[i] * inv_dist3);
        }
      }
    }
  }

  void IAS15::step(ParticleStore& particles) {
    size_t N = particles.N();
    if (N == 0) return;

    y0_.resize(N);
    v0_.resize(N);
    a0_.resize(N);
    for (size_t p = 0; p < N; ++p) {
      y0_[p] = particles.positions[p];
      v0_[p] = particles.velocities[p];
      a0_[p] = particles.accelerations[p];
    }

    a_sub_.resize(N * 8);
    for (size_t p = 0; p < N; ++p) {
      a_sub_[p * 8] = a0_[p];
    }

    g_.resize(N * 8);
    for (size_t p = 0; p < N; ++p) {
      g_[p * 8] = a0_[p];
    }

    for (size_t iter = 0; iter < max_iterations; ++iter) {
      int rr_idx = 0;
      for (size_t i = 1; i < 8; ++i) {
        for (size_t p = 0; p < N; ++p) {
          Vec3 sum_g{0, 0, 0};
          int rridx_local = rr_idx;
          for (size_t j = 1; j < i; ++j) {
            sum_g += g_[8 * p + j] * rr[rridx_local++];
          }
          g_[8 * p + i] = (a_sub_[p * 8 + i] - a0_[p] - sum_g) * (1.0 / h[i]);
        }
        rr_idx += i;      
      }

      for (size_t i = 1; i < 8; ++i) {
        double hn = h[i];
        double hdt = hn * dt_;
        double hdt2 = hdt * hdt;
        int c_idx = (i - 1) * i / 2;

        for (size_t p = 0; p < N; ++p) {
          Vec3 y_corr{0, 0, 0};
          Vec3 v_corr{0, 0, 0};

          int cid = c_idx;
          for (size_t j = 0; j < i; ++j) {
            double c_val = c[cid];
            double d_val = d[cid];
            y_corr += g_[p * 8 + j + 1] * c_val;
            v_corr += g_[p * 8 + j + 1] * d_val;
            cid++;
          }

          particles.positions[p] = y0_[p] + v0_[p] * hdt + a0_[p] * (0.5 * hdt2) + y_corr * (hdt2 * hdt);
          particles.velocities[p] = v0_[p] + a0_[p] * hdt + v_corr * hdt2;
        }

        std::vector<Vec3> current_accs(N);
        compute_accelerations(particles, current_accs);
        for (size_t p = 0; p < N; ++p) {
          a_sub_[p * 8 + i] = current_accs[p];
        }
      }
    }

    double max_delta = 0.0;
    for (size_t p = 0; p < N; ++p) {
      Vec3 sum_v{0, 0, 0};
      for (size_t i = 0; i < 8; ++i) {
        sum_v += a_sub_[p * 8 + i] * w[i];
      }
      particles.positions[p] = y0_[p] + v0_[p] * dt_ + a0_[p] * (0.5 * dt_ * dt_) + sum_v * (dt_ * dt_);
      particles.velocities[p] = v0_[p] + a0_[p] * dt_ + sum_v * dt_;
      particles.accelerations[p] = a_sub_[p * 8 + 7];

      double d_mag = g_[p * 8 + 7].mag();
      if (d_mag > max_delta) max_delta = d_mag;
    }

    if (max_delta > 0.0) {
      dt_ = dt_ * std::pow(epsilon_dtadjust / max_delta, 1.0 / 15.0) * 0.9;
    }
  }
} // namespace fnb
