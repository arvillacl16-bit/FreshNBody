// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VEC_HPP_
#define VEC_HPP_

#include <cmath>

#if defined(__GNUC__) || defined(__clang__)
#define restrict __restrict__
#elif defined(_MSC_VER)
#define restrict __restrict
#else
#define restrict
#endif

namespace fnb {
  struct Vec3;
  inline constexpr Vec3 operator-(Vec3 a, Vec3 b);

  struct Vec3 {
    double x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr explicit Vec3(double x_) : x(x_) {}
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    inline Vec3& operator+=(Vec3 other) {
      x += other.x;
      y += other.y;
      z += other.z;
      return *this;
    }

    inline Vec3& operator-=(Vec3 other) {
      x -= other.x;
      y -= other.y;
      z -= other.z;
      return *this;
    }

    inline Vec3& operator*=(double other) {
      x *= other;
      y *= other;
      z *= other;
      return *this;
    }

    inline Vec3& operator/=(double other) {
      double inv = 1. / other;
      operator*=(inv);
      return *this;
    }

    [[nodiscard]] inline constexpr double dot(Vec3 other) const { return x * other.x + y * other.y + z * other.z; }
    [[nodiscard]] inline constexpr double mag2() const { return dot(*this); }
    [[nodiscard]] inline double mag() const { return std::sqrt(mag2()); }
    [[nodiscard]] inline constexpr double dist2(Vec3 other) const { return (other - *this).mag2(); }
    [[nodiscard]] inline double dist(Vec3 other) const { return std::sqrt(dist2(other)); }
    [[nodiscard]] inline constexpr bool almost_eq(Vec3 other, double thres = 1e-9) const { return dist2(other) < thres * thres; }

    [[nodiscard]] inline constexpr Vec3 cross(Vec3 other) const {
      return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
    }

  private:
    template <typename T> inline constexpr T ce_abs(T a) { return a < 0 ? -a : a; }
  };

  [[nodiscard]] inline constexpr Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
  [[nodiscard]] inline constexpr Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
  [[nodiscard]] inline constexpr Vec3 operator*(Vec3 a, double b) { return {a.x * b, a.y * b, a.z * b}; }
  [[nodiscard]] inline constexpr Vec3 operator/(Vec3 a, double b) {
    double inv = 1. / b;
    return a * inv;
  }

  [[nodiscard]] inline constexpr Vec3 operator*(double a, Vec3 b) { return b * a; }
  [[nodiscard]] inline constexpr bool operator==(Vec3 a, Vec3 b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
  [[nodiscard]] inline constexpr bool operator!=(Vec3 a, Vec3 b) { return !(a == b); }
  [[nodiscard]] inline constexpr Vec3 operator+(Vec3 a) { return a; }
  [[nodiscard]] inline constexpr Vec3 operator-(Vec3 a) { return {-a.x, -a.y, -a.z}; }
} // namespace fnb
#pragma omp declare reduction(+ : fnb::Vec3 : omp_out += omp_in) initializer(omp_priv = fnb::Vec3{0, 0, 0})

#endif
