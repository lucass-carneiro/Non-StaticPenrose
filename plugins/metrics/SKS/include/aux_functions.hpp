#ifndef GRLENSING_SKS_AUX_FUNCTIONS
#define GRLENSING_SKS_AUX_FUNCTIONS

#include "../../../../include/tensor_types.hpp"

#include <cmath>
#include <cstddef>

namespace sks_aux {

enum class bh_tag : int { A = 1, B = 2 };

template <typename T> constexpr auto Power(T x, int n) {
  if (n == 0) {
    return T{1};
  } else if (n < 0) {
    return T{1} / Power(x, -n);
  } else {
    return x * Power(x, n - 1);
  }
}

constexpr auto Sqrt(auto x) {
  using std::sqrt;
  return sqrt(x);
}

inline auto r_KS(double a, double X, double Y, double Z) -> double {
  double part1 = Power(X, 2) + Power(Y, 2) + Power(Z, 2) - Power(a, 2);
  return Sqrt(part1 + Sqrt(4 * Power(a, 2) * Power(Z, 2) + Power(part1, 2))) / Sqrt(2);
}

inline auto d_r_KS_dX(double a, double X, double Y, double Z) -> double {
  return (X
          * Sqrt(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2)
                 + Sqrt(4 * Power(a, 2) * Power(Z, 2)
                        + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(4 * Power(a, 2) * Power(Z, 2)
                   + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2)));
}

inline auto d_r_KS_dY(double a, double X, double Y, double Z) -> double {
  return (Y
          * Sqrt(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2)
                 + Sqrt(4 * Power(a, 2) * Power(Z, 2)
                        + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(4 * Power(a, 2) * Power(Z, 2)
                   + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2)));
}

inline auto d_r_KS_dZ(double a, double X, double Y, double Z) -> double {
  return (Z
          * (1
             + (Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2))
                   / Sqrt(4 * Power(a, 2) * Power(Z, 2)
                          + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2))))
         / (Sqrt(2)
            * Sqrt(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2)
                   + Sqrt(4 * Power(a, 2) * Power(Z, 2)
                          + Power(-Power(a, 2) + Power(X, 2) + Power(Y, 2) + Power(Z, 2), 2))));
}

inline auto H_KS(double M, double a, double r, double Z) -> double {
  return (M * r) / (Power(r, 2) + Power(a, 2) * Power(Z / r, 2));
}

inline auto d_H_KS_dX(double M, double a, double dr_dX, double r, double Z) -> double {
  return -((M * Power(r, 2) * (-3 * Power(a, 2) * Power(Z, 2) + Power(r, 4)) * dr_dX)
           / Power(Power(a, 2) * Power(Z, 2) + Power(r, 4), 2));
}

inline auto d_H_KS_dY(double M, double a, double dr_dY, double r, double Z) -> double {
  return -((M * Power(r, 2) * (-3 * Power(a, 2) * Power(Z, 2) + Power(r, 4)) * dr_dY)
           / Power(Power(a, 2) * Power(Z, 2) + Power(r, 4), 2));
}

inline auto d_H_KS_dZ(double M, double a, double dr_dZ, double r, double Z) -> double {
  return (M * Power(r, 2)
          * (-2 * Power(a, 2) * Z * r + (3 * Power(a, 2) * Power(Z, 2) - Power(r, 4)) * dr_dZ))
         / Power(Power(a, 2) * Power(Z, 2) + Power(r, 4), 2);
}

inline auto l1_KS(double a, double r, double X, double Y) -> double {
  return (r * X + a * Y) / (Power(r, 2) + Power(a, 2));
}

inline auto d_l1_KS_dX(double a, double dr_dX, double r, double X, double Y) -> double {
  return (r * (Power(a, 2) + Power(r, 2)) + (Power(a, 2) * X - r * (2 * a * Y + X * r)) * dr_dX)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto d_l1_KS_dY(double a, double dr_dY, double r, double X, double Y) -> double {
  return (Power(a, 3) + a * Power(r, 2) + (Power(a, 2) * X - r * (2 * a * Y + X * r)) * dr_dY)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto d_l1_KS_dZ(double a, double dr_dZ, double r, double X, double Y) -> double {
  return ((Power(a, 2) * X - r * (2 * a * Y + X * r)) * dr_dZ)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto l2_KS(double a, double r, double X, double Y) -> double {
  return (r * Y - a * X) / (Power(r, 2) + Power(a, 2));
}

inline auto d_l2_KS_dX(double a, double dr_dX, double r, double X, double Y) -> double {
  return (-(a * (Power(a, 2) + Power(r, 2)))
          + (Power(a, 2) * Y + 2 * a * X * r - Y * Power(r, 2)) * dr_dX)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto d_l2_KS_dY(double a, double dr_dY, double r, double X, double Y) -> double {
  return (r * (Power(a, 2) + Power(r, 2))
          + (Power(a, 2) * Y + 2 * a * X * r - Y * Power(r, 2)) * dr_dY)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto d_l2_KS_dZ(double a, double dr_dZ, double r, double X, double Y) -> double {
  return ((Power(a, 2) * Y + 2 * a * X * r - Y * Power(r, 2)) * dr_dZ)
         / Power(Power(a, 2) + Power(r, 2), 2);
}

inline auto l3_KS(double r, double Z) -> double { return Z / r; }

inline auto d_l3_KS_dX(double r, double dr_dX, double Z) -> double {
  return -((Z * dr_dX) / Power(r, 2));
}

inline auto d_l3_KS_dY(double r, double dr_dY, double Z) -> double {
  return -((Z * dr_dY) / Power(r, 2));
}

inline auto d_l3_KS_dZ(double r, double dr_dZ, double Z) -> double {
  return (r - Z * dr_dZ) / Power(r, 2);
}

template <bh_tag tag> inline auto T(double b, double Omega, double t, double x, double y, double)
    -> double {
  using std::sqrt;

  return (8 * t * sqrt(4 - Power(b, 2) * Power(Omega, 2))
          - 2 * Power(-1, static_cast<int>(tag)) * b * y * Omega
                * sqrt(4 - Power(b, 2) * Power(Omega, 2)) * cos(3 * t * Omega)
          + 2 * Power(-1, static_cast<int>(tag)) * b * Omega * cos(t * Omega)
                * (-(y * sqrt(4 - Power(b, 2) * Power(Omega, 2))) + 4 * x * sin(2 * t * Omega))
          - 2 * b * Omega * sin(t * Omega)
                * (2 * sqrt(4 - Power(b, 2) * Power(Omega, 2))
                       * (Power(-1, static_cast<int>(tag)) * x + b * cos(t * Omega))
                       * cos(2 * t * Omega)
                   + 4 * Power(-1, static_cast<int>(tag)) * y * sin(2 * t * Omega))
          + 2 * Power(b, 2) * Omega * sin(4 * t * Omega))
         / 16;
}

template <bh_tag tag> inline auto X(double b, double Omega, double t, double x, double y, double)
    -> double {
  using std::sqrt;

  return (2 * x * sqrt(4 - Power(b, 2) * Power(Omega, 2))
          + Power(-1, static_cast<int>(tag)) * b * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                * cos(t * Omega)
          + 8 * x * Power(cos(t * Omega), 2)
          + 2 * Power(-1, static_cast<int>(tag)) * b * cos(3 * t * Omega)
          - 4 * y * sin(2 * t * Omega)
          - sqrt(4 - Power(b, 2) * Power(Omega, 2))
                * (2 * x * cos(2 * t * Omega)
                   + Power(-1, static_cast<int>(tag)) * b * cos(3 * t * Omega)
                   - 2 * y * sin(2 * t * Omega)))
         / 8;
}

template <bh_tag tag> inline auto Y(double b, double Omega, double t, double x, double y, double)
    -> double {
  using std::sqrt;
  return (y * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
          + sqrt(4 - Power(b, 2) * Power(Omega, 2))
                * (Power(-1, static_cast<int>(tag)) * b + 2 * x * cos(t * Omega)) * sin(t * Omega)
          + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)
                * (y + Power(-1, static_cast<int>(tag)) * b * sin(t * Omega))
          - 2 * x * sin(2 * t * Omega))
         / 4;
}

template <bh_tag tag> inline auto Z(double, double, double, double, double, double z) -> double {
  return z;
}

template <bh_tag tag> inline auto J(double b, double Omega, double t, double x, double y, double)
    -> grlensing::callable_matrix<4, 4, double> {

  using std::sqrt;

  grlensing::callable_matrix<4, 4, double> jac{};

  jac[0][0] = (4 * sqrt(4 - Power(b, 2) * Power(Omega, 2))
               + Power(-1, static_cast<int>(tag)) * b * x * Power(Omega, 2)
                     * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(t * Omega)
               - 3 * Power(-1, static_cast<int>(tag)) * b * x * Power(Omega, 2)
                     * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(3 * t * Omega)
               - b * Power(Omega, 2)
                     * (2 * b * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(4 * t * Omega)
                        - 2 * Power(-1, static_cast<int>(tag)) * y
                              * (2 * (-1 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                                 + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                                       * cos(2 * t * Omega))
                              * sin(t * Omega)))
              / 8.;

  jac[0][1] = -(Power(-1, static_cast<int>(tag)) * b * Omega
                * (-2 + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
                * sin(t * Omega))
              / 4.;

  jac[0][2] = -(Power(-1, static_cast<int>(tag)) * b * Omega * cos(t * Omega)
                * (2 + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)))
              / 4.;

  jac[1][0] = (Omega
                   * (Power(-1, static_cast<int>(tag)) * b
                          * (-4 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                      + 4 * x * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(t * Omega))
                   * sin(t * Omega)
               + Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)
                     * (2 * y + 3 * Power(-1, static_cast<int>(tag)) * b * sin(t * Omega)))
              / 4.;

  jac[1][1] = (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))
               - (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
              / 4.;

  jac[1][2] = ((-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 4.;

  jac[2][0] = (Omega
               * (Power(-1, static_cast<int>(tag)) * b
                      * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(t * Omega)
                  + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                        * (4 * x * cos(2 * t * Omega)
                           + 3 * Power(-1, static_cast<int>(tag)) * b * cos(3 * t * Omega)
                           - 4 * y * sin(2 * t * Omega))))
              / 8.;

  jac[2][1] = ((-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 4.;

  jac[2][2] = (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))
               + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
              / 4.;

  jac[3][3] = 1;

  return jac;
}

template <bh_tag tag> inline auto dJ(double b, double Omega, double t, double x, double y, double)
    -> grlensing::callable_3tensor<4, 4, 4, double> {

  using std::sqrt;

  grlensing::callable_3tensor<4, 4, 4, double> djac{};

  djac[0][0][0]
      = (b * Power(Omega, 3)
         * (Power(-1, static_cast<int>(tag)) * y * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                * cos(t * Omega)
            + 9 * Power(-1, static_cast<int>(tag)) * y
                  * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(3 * t * Omega)
            + Power(-1, 1 + static_cast<int>(tag)) * x * sqrt(4 - Power(b, 2) * Power(Omega, 2))
                  * sin(t * Omega)
            + 9 * Power(-1, static_cast<int>(tag)) * x * sqrt(4 - Power(b, 2) * Power(Omega, 2))
                  * sin(3 * t * Omega)
            - 2 * Power(-1, static_cast<int>(tag)) * x * (sin(t * Omega) + 9 * sin(3 * t * Omega))
            - 16 * b * sin(4 * t * Omega)
            + 8 * b * sqrt(4 - Power(b, 2) * Power(Omega, 2)) * sin(4 * t * Omega)))
        / 8.;

  djac[0][0][1] = -(Power(-1, static_cast<int>(tag)) * b * Power(Omega, 2) * cos(t * Omega)
                    * (2 - 2 * sqrt(4 - Power(b, 2) * Power(Omega, 2))
                       + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)))
                  / 4.;

  djac[0][0][2] = (Power(-1, static_cast<int>(tag)) * b * Power(Omega, 2)
                   * (2 * (-1 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                      + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
                   * sin(t * Omega))
                  / 4.;

  djac[0][1][0] = -(Power(-1, static_cast<int>(tag)) * b * Power(Omega, 2) * cos(t * Omega)
                    * (2 - 2 * sqrt(4 - Power(b, 2) * Power(Omega, 2))
                       + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)))
                  / 4.;

  djac[0][2][0] = (Power(-1, static_cast<int>(tag)) * b * Power(Omega, 2)
                   * (2 * (-1 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                      + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
                   * sin(t * Omega))
                  / 4.;

  djac[1][0][0] = -(Power(Omega, 2)
                    * (Power(-1, static_cast<int>(tag)) * b
                           * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(t * Omega)
                       - (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                             * (8 * x * cos(2 * t * Omega)
                                + 9 * Power(-1, static_cast<int>(tag)) * b * cos(3 * t * Omega)
                                - 8 * y * sin(2 * t * Omega))))
                  / 8.;

  djac[1][0][1]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 2.;

  djac[1][0][2]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)) / 2.;

  djac[1][1][0]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 2.;

  djac[1][2][0]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)) / 2.;

  djac[2][0][0] = -(Power(Omega, 2)
                    * (8 * y * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)
                       + Power(-1, static_cast<int>(tag)) * b
                             * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(t * Omega)
                       + (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                             * (8 * x * sin(2 * t * Omega)
                                + 9 * Power(-1, static_cast<int>(tag)) * b * sin(3 * t * Omega))))
                  / 8.;

  djac[2][0][1]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)) / 2.;

  djac[2][0][2]
      = -(Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 2.;

  djac[2][1][0]
      = (Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega)) / 2.;

  djac[2][2][0]
      = -(Omega * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * sin(2 * t * Omega)) / 2.;

  return djac;
}

inline auto inverse_symmetric_4x4(const grlensing::callable_matrix<4, 4, double> &matrix)
    -> grlensing::callable_matrix<4, 4, double> {
  grlensing::callable_matrix<4, 4, double> inverse{};

  const auto det
      = -(matrix(0, 0) * Power(matrix(1, 3), 2) * matrix(2, 2))
        + Power(matrix(0, 3), 2) * (Power(matrix(1, 2), 2) - matrix(1, 1) * matrix(2, 2))
        + 2 * matrix(0, 0) * matrix(1, 2) * matrix(1, 3) * matrix(2, 3)
        + Power(matrix(0, 1), 2) * Power(matrix(2, 3), 2)
        - matrix(0, 0) * matrix(1, 1) * Power(matrix(2, 3), 2)
        + 2 * matrix(0, 3)
              * (matrix(0, 2) * (-(matrix(1, 2) * matrix(1, 3)) + matrix(1, 1) * matrix(2, 3))
                 + matrix(0, 1) * (matrix(1, 3) * matrix(2, 2) - matrix(1, 2) * matrix(2, 3)))
        - (matrix(0, 0) * Power(matrix(1, 2), 2)
           + (Power(matrix(0, 1), 2) - matrix(0, 0) * matrix(1, 1)) * matrix(2, 2))
              * matrix(3, 3)
        + Power(matrix(0, 2), 2) * (Power(matrix(1, 3), 2) - matrix(1, 1) * matrix(3, 3))
        + 2 * matrix(0, 1) * matrix(0, 2)
              * (-(matrix(1, 3) * matrix(2, 3)) + matrix(1, 2) * matrix(3, 3));

  inverse[0][0]
      = (-(Power(matrix(1, 3), 2) * matrix(2, 2)) + 2 * matrix(1, 2) * matrix(1, 3) * matrix(2, 3)
         - Power(matrix(1, 2), 2) * matrix(3, 3)
         + matrix(1, 1) * (-Power(matrix(2, 3), 2) + matrix(2, 2) * matrix(3, 3)))
        / det;

  inverse[0][1] = (matrix(0, 3) * (matrix(1, 3) * matrix(2, 2) - matrix(1, 2) * matrix(2, 3))
                   + matrix(0, 2) * (-(matrix(1, 3) * matrix(2, 3)) + matrix(1, 2) * matrix(3, 3))
                   + matrix(0, 1) * (Power(matrix(2, 3), 2) - matrix(2, 2) * matrix(3, 3)))
                  / det;

  inverse[0][2] = (matrix(0, 3) * (-(matrix(1, 2) * matrix(1, 3)) + matrix(1, 1) * matrix(2, 3))
                   + matrix(0, 2) * (Power(matrix(1, 3), 2) - matrix(1, 1) * matrix(3, 3))
                   + matrix(0, 1) * (-(matrix(1, 3) * matrix(2, 3)) + matrix(1, 2) * matrix(3, 3)))
                  / det;

  inverse[0][3] = (matrix(0, 3) * (Power(matrix(1, 2), 2) - matrix(1, 1) * matrix(2, 2))
                   + matrix(0, 2) * (-(matrix(1, 2) * matrix(1, 3)) + matrix(1, 1) * matrix(2, 3))
                   + matrix(0, 1) * (matrix(1, 3) * matrix(2, 2) - matrix(1, 2) * matrix(2, 3)))
                  / det;

  inverse[1][1]
      = (-(Power(matrix(0, 3), 2) * matrix(2, 2)) + 2 * matrix(0, 2) * matrix(0, 3) * matrix(2, 3)
         - Power(matrix(0, 2), 2) * matrix(3, 3)
         + matrix(0, 0) * (-Power(matrix(2, 3), 2) + matrix(2, 2) * matrix(3, 3)))
        / det;

  inverse[1][2] = (Power(matrix(0, 3), 2) * matrix(1, 2)
                   - matrix(0, 3) * (matrix(0, 2) * matrix(1, 3) + matrix(0, 1) * matrix(2, 3))
                   + matrix(0, 1) * matrix(0, 2) * matrix(3, 3)
                   + matrix(0, 0) * (matrix(1, 3) * matrix(2, 3) - matrix(1, 2) * matrix(3, 3)))
                  / det;

  inverse[1][3]
      = (Power(matrix(0, 2), 2) * matrix(1, 3) + matrix(0, 1) * matrix(0, 3) * matrix(2, 2)
         - matrix(0, 2) * (matrix(0, 3) * matrix(1, 2) + matrix(0, 1) * matrix(2, 3))
         + matrix(0, 0) * (-(matrix(1, 3) * matrix(2, 2)) + matrix(1, 2) * matrix(2, 3)))
        / det;

  inverse[2][2]
      = (-(Power(matrix(0, 3), 2) * matrix(1, 1)) + 2 * matrix(0, 1) * matrix(0, 3) * matrix(1, 3)
         - Power(matrix(0, 1), 2) * matrix(3, 3)
         + matrix(0, 0) * (-Power(matrix(1, 3), 2) + matrix(1, 1) * matrix(3, 3)))
        / det;

  inverse[2][3] = (matrix(1, 2) * (-(matrix(0, 1) * matrix(0, 3)) + matrix(0, 0) * matrix(1, 3))
                   + matrix(0, 2) * (matrix(0, 3) * matrix(1, 1) - matrix(0, 1) * matrix(1, 3))
                   + (Power(matrix(0, 1), 2) - matrix(0, 0) * matrix(1, 1)) * matrix(2, 3))
                  / det;

  inverse[3][3]
      = (-(Power(matrix(0, 2), 2) * matrix(1, 1)) + 2 * matrix(0, 1) * matrix(0, 2) * matrix(1, 2)
         - Power(matrix(0, 1), 2) * matrix(2, 2)
         + matrix(0, 0) * (-Power(matrix(1, 2), 2) + matrix(1, 1) * matrix(2, 2)))
        / det;

  inverse[1][0] = inverse[0][1];
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];
  inverse[3][0] = inverse[0][3];
  inverse[3][1] = inverse[1][3];
  inverse[3][2] = inverse[2][3];

  return inverse;
}

inline auto inverse_symmetric_3x3(const grlensing::callable_matrix<3, 3, double> &matrix)
    -> grlensing::callable_matrix<3, 3, double> {
  grlensing::callable_matrix<3, 3, double> inverse{};

  const auto det = -(Power(matrix(0, 2), 2) * matrix(1, 1))
                   + 2 * matrix(0, 1) * matrix(0, 2) * matrix(1, 2)
                   - Power(matrix(0, 1), 2) * matrix(2, 2)
                   + matrix(0, 0) * (-Power(matrix(1, 2), 2) + matrix(1, 1) * matrix(2, 2));

  inverse[0][0] = (-Power(matrix(1, 2), 2) + matrix(1, 1) * matrix(2, 2)) / det;

  inverse[0][1] = (matrix(0, 2) * matrix(1, 2) - matrix(0, 1) * matrix(2, 2)) / det;

  inverse[0][2] = (-(matrix(0, 2) * matrix(1, 1)) + matrix(0, 1) * matrix(1, 2)) / det;

  inverse[1][1] = (-Power(matrix(0, 2), 2) + matrix(0, 0) * matrix(2, 2)) / det;

  inverse[1][2] = (matrix(0, 1) * matrix(0, 2) - matrix(0, 0) * matrix(1, 2)) / det;

  inverse[2][2] = (-Power(matrix(0, 1), 2) + matrix(0, 0) * matrix(1, 1)) / det;

  inverse[1][0] = inverse[0][1];
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];

  return inverse;
}

inline auto llg_SKS(double M1, double M2, double a1, double a2, double b, double t, double x,
                    double y, double z) -> grlensing::callable_matrix<4, 4, double> {

  double Omega = sqrt((M1 + M2) / Power(b, 3));

  const auto X1 = X<bh_tag::A>(b, Omega, t, x, y, z);
  const auto Y1 = Y<bh_tag::A>(b, Omega, t, x, y, z);
  const auto Z1 = Z<bh_tag::A>(b, Omega, t, x, y, z);

  const auto X2 = X<bh_tag::B>(b, Omega, t, x, y, z);
  const auto Y2 = Y<bh_tag::B>(b, Omega, t, x, y, z);
  const auto Z2 = Z<bh_tag::B>(b, Omega, t, x, y, z);

  const auto rKS1 = r_KS(a1, X1, Y1, Z1);
  const auto rKS2 = r_KS(a2, X2, Y2, Z2);

  const auto H1 = H_KS(M1, a1, rKS1, Z1);
  const auto H2 = H_KS(M2, a2, rKS2, Z2);

  const auto l1_1 = l1_KS(a1, rKS1, X1, Y1);
  const auto l1_2 = l1_KS(a2, rKS2, X2, Y2);

  const auto l2_1 = l2_KS(a1, rKS1, X1, Y1);
  const auto l2_2 = l2_KS(a2, rKS2, X2, Y2);

  const auto l3_1 = l3_KS(rKS1, Z1);
  const auto l3_2 = l3_KS(rKS2, Z2);

  const auto jac1 = J<bh_tag::A>(b, Omega, t, x, y, z);
  const auto jac2 = J<bh_tag::B>(b, Omega, t, x, y, z);

  grlensing::callable_matrix<4, 4, double> llg{};

  llg[0][0]
      = -1
        + 2 * H1 * Power(jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1, 2)
        + 2 * H2 * Power(jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2, 2);

  llg[0][1] = 2 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                  * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
              + 2 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                    * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2);

  llg[0][2] = 2 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                  * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
              + 2 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                    * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2);

  llg[0][3] = 2 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                  * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
              + 2 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                    * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2);

  llg[1][1]
      = 1
        + 2 * H1 * Power(jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1, 2)
        + 2 * H2 * Power(jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2, 2);

  llg[1][2] = 2 * H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                  * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
              + 2 * H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                    * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2);

  llg[1][3] = 2 * H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                  * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
              + 2 * H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                    * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2);

  llg[2][2]
      = 1
        + 2 * H1 * Power(jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1, 2)
        + 2 * H2 * Power(jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2, 2);

  llg[2][3] = 2 * H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                  * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
              + 2 * H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                    * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2);

  llg[3][3]
      = 1
        + 2 * H1 * Power(jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1, 2)
        + 2 * H2 * Power(jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2, 2);

  llg[1][0] = llg[0][1];
  llg[2][0] = llg[0][2];
  llg[2][1] = llg[1][2];
  llg[3][0] = llg[0][3];
  llg[3][1] = llg[1][3];
  llg[3][2] = llg[2][3];

  return llg;
}

inline auto dllg_SKS(double M1, double M2, double a1, double a2, double b, double t, double x,
                     double y, double z) -> grlensing::callable_3tensor<4, 4, 4, double> {

  double Omega = sqrt((M1 + M2) / Power(b, 3));

  const auto X1 = X<bh_tag::A>(b, Omega, t, x, y, z);
  const auto Y1 = Y<bh_tag::A>(b, Omega, t, x, y, z);
  const auto Z1 = Z<bh_tag::A>(b, Omega, t, x, y, z);

  const auto X2 = X<bh_tag::B>(b, Omega, t, x, y, z);
  const auto Y2 = Y<bh_tag::B>(b, Omega, t, x, y, z);
  const auto Z2 = Z<bh_tag::B>(b, Omega, t, x, y, z);

  const auto rKS1 = r_KS(a1, X1, Y1, Z1);
  const auto rKS2 = r_KS(a2, X2, Y2, Z2);

  const auto H1 = H_KS(M1, a1, rKS1, Z1);
  const auto H2 = H_KS(M2, a2, rKS2, Z2);

  const auto l1_1 = l1_KS(a1, rKS1, X1, Y1);
  const auto l1_2 = l1_KS(a2, rKS2, X2, Y2);

  const auto l2_1 = l2_KS(a1, rKS1, X1, Y1);
  const auto l2_2 = l2_KS(a2, rKS2, X2, Y2);

  const auto l3_1 = l3_KS(rKS1, Z1);
  const auto l3_2 = l3_KS(rKS2, Z2);

  const auto jac1 = J<bh_tag::A>(b, Omega, t, x, y, z);
  const auto jac2 = J<bh_tag::B>(b, Omega, t, x, y, z);

  const auto drKS1_dX1 = d_r_KS_dX(a1, X1, Y1, Z1);
  const auto drKS1_dY1 = d_r_KS_dY(a1, X1, Y1, Z1);
  const auto drKS1_dZ1 = d_r_KS_dZ(a1, X1, Y1, Z1);

  const auto drKS2_dX2 = d_r_KS_dX(a2, X2, Y2, Z2);
  const auto drKS2_dY2 = d_r_KS_dY(a2, X2, Y2, Z2);
  const auto drKS2_dZ2 = d_r_KS_dZ(a2, X2, Y2, Z2);

  const auto dH1_dX1 = d_H_KS_dX(M1, a1, drKS1_dX1, rKS1, Z1);
  const auto dH1_dY1 = d_H_KS_dY(M1, a1, drKS1_dY1, rKS1, Z1);
  const auto dH1_dZ1 = d_H_KS_dZ(M1, a1, drKS1_dZ1, rKS1, Z1);

  const auto dH2_dX2 = d_H_KS_dX(M2, a2, drKS2_dX2, rKS2, Z2);
  const auto dH2_dY2 = d_H_KS_dY(M2, a2, drKS2_dY2, rKS2, Z2);
  const auto dH2_dZ2 = d_H_KS_dZ(M2, a2, drKS2_dZ2, rKS2, Z2);

  const auto dl1_1_dX1 = d_l1_KS_dX(a1, drKS1_dX1, rKS1, X1, Y1);
  const auto dl1_1_dY1 = d_l1_KS_dY(a1, drKS1_dY1, rKS1, X1, Y1);
  const auto dl1_1_dZ1 = d_l1_KS_dZ(a1, drKS1_dZ1, rKS1, X1, Y1);

  const auto dl2_1_dX1 = d_l2_KS_dX(a1, drKS1_dX1, rKS1, X1, Y1);
  const auto dl2_1_dY1 = d_l2_KS_dY(a1, drKS1_dY1, rKS1, X1, Y1);
  const auto dl2_1_dZ1 = d_l2_KS_dZ(a1, drKS1_dZ1, rKS1, X1, Y1);

  const auto dl3_1_dX1 = d_l3_KS_dX(rKS1, drKS1_dX1, Z1);
  const auto dl3_1_dY1 = d_l3_KS_dY(rKS1, drKS1_dY1, Z1);
  const auto dl3_1_dZ1 = d_l3_KS_dZ(rKS1, drKS1_dZ1, Z1);

  const auto dl1_2_dX2 = d_l1_KS_dX(a2, drKS2_dX2, rKS2, X2, Y2);
  const auto dl1_2_dY2 = d_l1_KS_dY(a2, drKS2_dY2, rKS2, X2, Y2);
  const auto dl1_2_dZ2 = d_l1_KS_dZ(a2, drKS2_dZ2, rKS2, X2, Y2);

  const auto dl2_2_dX2 = d_l2_KS_dX(a2, drKS2_dX2, rKS2, X2, Y2);
  const auto dl2_2_dY2 = d_l2_KS_dY(a2, drKS2_dY2, rKS2, X2, Y2);
  const auto dl2_2_dZ2 = d_l2_KS_dZ(a2, drKS2_dZ2, rKS2, X2, Y2);

  const auto dl3_2_dX2 = d_l3_KS_dX(rKS2, drKS2_dX2, Z2);
  const auto dl3_2_dY2 = d_l3_KS_dY(rKS2, drKS2_dY2, Z2);
  const auto dl3_2_dZ2 = d_l3_KS_dZ(rKS2, drKS2_dZ2, Z2);

  const auto djac1 = dJ<bh_tag::A>(b, Omega, t, x, y, z);
  const auto djac2 = dJ<bh_tag::B>(b, Omega, t, x, y, z);

  grlensing::callable_3tensor<4, 4, 4, double> dllg{};

  dllg[0][0][0]
      = 4 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
            * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
               + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
               + Power(jac1(2, 0), 2) * dl2_1_dY1
               + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
               + Power(jac1(1, 0), 2) * dl1_1_dX1
               + jac1(1, 0)
                     * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                        + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
        + 4 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
              * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                 + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                 + Power(jac2(2, 0), 2) * dl2_2_dY2
                 + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                 + Power(jac2(1, 0), 2) * dl1_2_dX2
                 + jac2(1, 0)
                       * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                          + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
        + 2 * Power(jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1, 2)
              * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
        + 2 * Power(jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2, 2)
              * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2);
  dllg[1][0][0]
      = 4 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
            * (djac1(0, 0, 1) + djac1(1, 0, 1) * l1_1 + djac1(2, 0, 1) * l2_1
               + djac1(3, 0, 1) * l3_1
               + jac1(1, 0)
                     * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
               + jac1(2, 0)
                     * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
               + jac1(3, 0)
                     * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1 + jac1(1, 1) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
              * (djac2(0, 0, 1) + djac2(1, 0, 1) * l1_2 + djac2(2, 0, 1) * l2_2
                 + djac2(3, 0, 1) * l3_2
                 + jac2(1, 0)
                       * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2 + jac2(1, 1) * dl1_2_dX2)
                 + jac2(2, 0)
                       * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2 + jac2(1, 1) * dl2_2_dX2)
                 + jac2(3, 0)
                       * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2 + jac2(1, 1) * dl3_2_dX2))
        + 2 * Power(jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1, 2)
              * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
        + 2 * Power(jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2, 2)
              * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2);
  dllg[2][0][0]
      = 4 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
            * (djac1(0, 0, 2) + djac1(1, 0, 2) * l1_1 + djac1(2, 0, 2) * l2_1
               + djac1(3, 0, 2) * l3_1
               + jac1(1, 0)
                     * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
               + jac1(2, 0)
                     * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
               + jac1(3, 0)
                     * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1 + jac1(1, 2) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
              * (djac2(0, 0, 2) + djac2(1, 0, 2) * l1_2 + djac2(2, 0, 2) * l2_2
                 + djac2(3, 0, 2) * l3_2
                 + jac2(1, 0)
                       * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2 + jac2(1, 2) * dl1_2_dX2)
                 + jac2(2, 0)
                       * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2 + jac2(1, 2) * dl2_2_dX2)
                 + jac2(3, 0)
                       * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2 + jac2(1, 2) * dl3_2_dX2))
        + 2 * Power(jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1, 2)
              * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
        + 2 * Power(jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2, 2)
              * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2);
  dllg[3][0][0]
      = 4 * H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
            * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
               + djac1(3, 0, 3) * l3_1
               + jac1(1, 0)
                     * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
               + jac1(2, 0)
                     * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
               + jac1(3, 0)
                     * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1 + jac1(1, 3) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
              * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                 + djac2(3, 0, 3) * l3_2
                 + jac2(1, 0)
                       * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2 + jac2(1, 3) * dl1_2_dX2)
                 + jac2(2, 0)
                       * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2 + jac2(1, 3) * dl2_2_dX2)
                 + jac2(3, 0)
                       * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2 + jac2(1, 3) * dl3_2_dX2))
        + 2 * Power(jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1, 2)
              * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
        + 2 * Power(jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2, 2)
              * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2);

  dllg[0][0][1]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 1) + djac1(1, 0, 1) * l1_1 + djac1(2, 0, 1) * l2_1
                  + djac1(3, 0, 1) * l3_1
                  + jac1(1, 1)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 1)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 1)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 1) + djac2(1, 0, 1) * l1_2 + djac2(2, 0, 1) * l2_2
                    + djac2(3, 0, 1) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][0][1]
      = 2
        * (H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
               * (djac1(0, 0, 1) + djac1(1, 0, 1) * l1_1 + djac1(2, 0, 1) * l2_1
                  + djac1(3, 0, 1) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 1, 1) + djac1(1, 1, 1) * l1_1 + djac1(2, 1, 1) * l2_1
                    + djac1(3, 1, 1) * l3_1 + Power(jac1(3, 1), 2) * dl3_1_dZ1
                    + Power(jac1(2, 1), 2) * dl2_1_dY1
                    + jac1(2, 1) * jac1(3, 1) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 1), 2) * dl1_1_dX1
                    + jac1(1, 1)
                          * (jac1(2, 1) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 1) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 0, 1) + djac2(1, 0, 1) * l1_2 + djac2(2, 0, 1) * l2_2
                    + djac2(3, 0, 1) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 1, 1) + djac2(1, 1, 1) * l1_2 + djac2(2, 1, 1) * l2_2
                    + djac2(3, 1, 1) * l3_2 + Power(jac2(3, 1), 2) * dl3_2_dZ2
                    + Power(jac2(2, 1), 2) * dl2_2_dY2
                    + jac2(2, 1) * jac2(3, 1) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 1), 2) * dl1_2_dX2
                    + jac2(1, 1)
                          * (jac2(2, 1) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 1) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][0][1]
      = 2
        * (H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
               * (djac1(0, 0, 2) + djac1(1, 0, 2) * l1_1 + djac1(2, 0, 2) * l2_1
                  + djac1(3, 0, 2) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                    + djac1(3, 1, 2) * l3_1
                    + jac1(1, 1)
                          * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1
                             + jac1(1, 2) * dl1_1_dX1)
                    + jac1(2, 1)
                          * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1
                             + jac1(1, 2) * dl2_1_dX1)
                    + jac1(3, 1)
                          * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                             + jac1(1, 2) * dl3_1_dX1))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 0, 2) + djac2(1, 0, 2) * l1_2 + djac2(2, 0, 2) * l2_2
                    + djac2(3, 0, 2) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][0][1]
      = 2
        * (H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                    + djac1(3, 1, 3) * l3_1
                    + jac1(1, 1)
                          * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1
                             + jac1(1, 3) * dl1_1_dX1)
                    + jac1(2, 1)
                          * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1
                             + jac1(1, 3) * dl2_1_dX1)
                    + jac1(3, 1)
                          * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                             + jac1(1, 3) * dl3_1_dX1))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][0][2]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 2) + djac1(1, 0, 2) * l1_1 + djac1(2, 0, 2) * l2_1
                  + djac1(3, 0, 2) * l3_1
                  + jac1(1, 2)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 2)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 2)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 2) + djac2(1, 0, 2) * l1_2 + djac2(2, 0, 2) * l2_2
                    + djac2(3, 0, 2) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][0][2]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 0, 1) + djac1(1, 0, 1) * l1_1 + djac1(2, 0, 1) * l2_1
                  + djac1(3, 0, 1) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                    + djac1(3, 1, 2) * l3_1
                    + jac1(1, 2)
                          * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1
                             + jac1(1, 1) * dl1_1_dX1)
                    + jac1(2, 2)
                          * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1
                             + jac1(1, 1) * dl2_1_dX1)
                    + jac1(3, 2)
                          * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                             + jac1(1, 1) * dl3_1_dX1))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 0, 1) + djac2(1, 0, 1) * l1_2 + djac2(2, 0, 1) * l2_2
                    + djac2(3, 0, 1) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][0][2]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 0, 2) + djac1(1, 0, 2) * l1_1 + djac1(2, 0, 2) * l2_1
                  + djac1(3, 0, 2) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 2, 2) + djac1(1, 2, 2) * l1_1 + djac1(2, 2, 2) * l2_1
                    + djac1(3, 2, 2) * l3_1 + Power(jac1(3, 2), 2) * dl3_1_dZ1
                    + Power(jac1(2, 2), 2) * dl2_1_dY1
                    + jac1(2, 2) * jac1(3, 2) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 2), 2) * dl1_1_dX1
                    + jac1(1, 2)
                          * (jac1(2, 2) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 2) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 0, 2) + djac2(1, 0, 2) * l1_2 + djac2(2, 0, 2) * l2_2
                    + djac2(3, 0, 2) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 2, 2) + djac2(1, 2, 2) * l1_2 + djac2(2, 2, 2) * l2_2
                    + djac2(3, 2, 2) * l3_2 + Power(jac2(3, 2), 2) * dl3_2_dZ2
                    + Power(jac2(2, 2), 2) * dl2_2_dY2
                    + jac2(2, 2) * jac2(3, 2) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 2), 2) * dl1_2_dX2
                    + jac2(1, 2)
                          * (jac2(2, 2) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 2) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][0][2]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                    + djac1(3, 2, 3) * l3_1
                    + jac1(1, 2)
                          * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1
                             + jac1(1, 3) * dl1_1_dX1)
                    + jac1(2, 2)
                          * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1
                             + jac1(1, 3) * dl2_1_dX1)
                    + jac1(3, 2)
                          * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                             + jac1(1, 3) * dl3_1_dX1))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][0][3]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][0][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 0, 1) + djac1(1, 0, 1) * l1_1 + djac1(2, 0, 1) * l2_1
                  + djac1(3, 0, 1) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                    + djac1(3, 1, 3) * l3_1
                    + jac1(1, 3)
                          * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1
                             + jac1(1, 1) * dl1_1_dX1)
                    + jac1(2, 3)
                          * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1
                             + jac1(1, 1) * dl2_1_dX1)
                    + jac1(3, 3)
                          * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                             + jac1(1, 1) * dl3_1_dX1))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 1) + djac2(1, 0, 1) * l1_2 + djac2(2, 0, 1) * l2_2
                    + djac2(3, 0, 1) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][0][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 0, 2) + djac1(1, 0, 2) * l1_1 + djac1(2, 0, 2) * l2_1
                  + djac1(3, 0, 2) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                    + djac1(3, 2, 3) * l3_1
                    + jac1(1, 3)
                          * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1
                             + jac1(1, 2) * dl1_1_dX1)
                    + jac1(2, 3)
                          * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1
                             + jac1(1, 2) * dl2_1_dX1)
                    + jac1(3, 3)
                          * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                             + jac1(1, 2) * dl3_1_dX1))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 2) + djac2(1, 0, 2) * l1_2 + djac2(2, 0, 2) * l2_2
                    + djac2(3, 0, 2) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][0][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 0)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 0)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 0)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (djac1(0, 3, 3) + djac1(1, 3, 3) * l1_1 + djac1(2, 3, 3) * l2_1
                    + djac1(3, 3, 3) * l3_1 + Power(jac1(3, 3), 2) * dl3_1_dZ1
                    + Power(jac1(2, 3), 2) * dl2_1_dY1
                    + jac1(2, 3) * jac1(3, 3) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 3), 2) * dl1_1_dX1
                    + jac1(1, 3)
                          * (jac1(2, 3) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 3) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 0)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 0)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 0)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 3, 3) + djac2(1, 3, 3) * l1_2 + djac2(2, 3, 3) * l2_2
                    + djac2(3, 3, 3) * l3_2 + Power(jac2(3, 3), 2) * dl3_2_dZ2
                    + Power(jac2(2, 3), 2) * dl2_2_dY2
                    + jac2(2, 3) * jac2(3, 3) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 3), 2) * dl1_2_dX2
                    + jac2(1, 3)
                          * (jac2(2, 3) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 3) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][1][1]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][1][1]
      = 4 * H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
            * (djac1(0, 1, 1) + djac1(1, 1, 1) * l1_1 + djac1(2, 1, 1) * l2_1
               + djac1(3, 1, 1) * l3_1 + Power(jac1(3, 1), 2) * dl3_1_dZ1
               + Power(jac1(2, 1), 2) * dl2_1_dY1
               + jac1(2, 1) * jac1(3, 1) * (dl2_1_dZ1 + dl3_1_dY1)
               + Power(jac1(1, 1), 2) * dl1_1_dX1
               + jac1(1, 1)
                     * (jac1(2, 1) * (dl1_1_dY1 + dl2_1_dX1)
                        + jac1(3, 1) * (dl1_1_dZ1 + dl3_1_dX1)))
        + 4 * H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
              * (djac2(0, 1, 1) + djac2(1, 1, 1) * l1_2 + djac2(2, 1, 1) * l2_2
                 + djac2(3, 1, 1) * l3_2 + Power(jac2(3, 1), 2) * dl3_2_dZ2
                 + Power(jac2(2, 1), 2) * dl2_2_dY2
                 + jac2(2, 1) * jac2(3, 1) * (dl2_2_dZ2 + dl3_2_dY2)
                 + Power(jac2(1, 1), 2) * dl1_2_dX2
                 + jac2(1, 1)
                       * (jac2(2, 1) * (dl1_2_dY2 + dl2_2_dX2)
                          + jac2(3, 1) * (dl1_2_dZ2 + dl3_2_dX2)))
        + 2 * Power(jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1, 2)
              * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
        + 2 * Power(jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2, 2)
              * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2);
  dllg[2][1][1]
      = 4 * H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
            * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
               + djac1(3, 1, 2) * l3_1
               + jac1(1, 1)
                     * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
               + jac1(2, 1)
                     * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
               + jac1(3, 1)
                     * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1 + jac1(1, 2) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
              * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                 + djac2(3, 1, 2) * l3_2
                 + jac2(1, 1)
                       * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2 + jac2(1, 2) * dl1_2_dX2)
                 + jac2(2, 1)
                       * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2 + jac2(1, 2) * dl2_2_dX2)
                 + jac2(3, 1)
                       * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2 + jac2(1, 2) * dl3_2_dX2))
        + 2 * Power(jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1, 2)
              * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
        + 2 * Power(jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2, 2)
              * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2);
  dllg[3][1][1]
      = 4 * H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
            * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
               + djac1(3, 1, 3) * l3_1
               + jac1(1, 1)
                     * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
               + jac1(2, 1)
                     * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
               + jac1(3, 1)
                     * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1 + jac1(1, 3) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
              * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                 + djac2(3, 1, 3) * l3_2
                 + jac2(1, 1)
                       * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2 + jac2(1, 3) * dl1_2_dX2)
                 + jac2(2, 1)
                       * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2 + jac2(1, 3) * dl2_2_dX2)
                 + jac2(3, 1)
                       * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2 + jac2(1, 3) * dl3_2_dX2))
        + 2 * Power(jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1, 2)
              * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
        + 2 * Power(jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2, 2)
              * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2);

  dllg[0][1][2]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][1][2]
      = 2
        * (H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
               * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                  + djac1(3, 1, 2) * l3_1
                  + jac1(1, 2)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 2)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 2)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (djac1(0, 1, 1) + djac1(1, 1, 1) * l1_1 + djac1(2, 1, 1) * l2_1
                    + djac1(3, 1, 1) * l3_1 + Power(jac1(3, 1), 2) * dl3_1_dZ1
                    + Power(jac1(2, 1), 2) * dl2_1_dY1
                    + jac1(2, 1) * jac1(3, 1) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 1), 2) * dl1_1_dX1
                    + jac1(1, 1)
                          * (jac1(2, 1) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 1) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 1, 1) + djac2(1, 1, 1) * l1_2 + djac2(2, 1, 1) * l2_2
                    + djac2(3, 1, 1) * l3_2 + Power(jac2(3, 1), 2) * dl3_2_dZ2
                    + Power(jac2(2, 1), 2) * dl2_2_dY2
                    + jac2(2, 1) * jac2(3, 1) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 1), 2) * dl1_2_dX2
                    + jac2(1, 1)
                          * (jac2(2, 1) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 1) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][1][2]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                  + djac1(3, 1, 2) * l3_1
                  + jac1(1, 1)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 1)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 1)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (djac1(0, 2, 2) + djac1(1, 2, 2) * l1_1 + djac1(2, 2, 2) * l2_1
                    + djac1(3, 2, 2) * l3_1 + Power(jac1(3, 2), 2) * dl3_1_dZ1
                    + Power(jac1(2, 2), 2) * dl2_1_dY1
                    + jac1(2, 2) * jac1(3, 2) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 2), 2) * dl1_1_dX1
                    + jac1(1, 2)
                          * (jac1(2, 2) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 2) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 2, 2) + djac2(1, 2, 2) * l1_2 + djac2(2, 2, 2) * l2_2
                    + djac2(3, 2, 2) * l3_2 + Power(jac2(3, 2), 2) * dl3_2_dZ2
                    + Power(jac2(2, 2), 2) * dl2_2_dY2
                    + jac2(2, 2) * jac2(3, 2) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 2), 2) * dl1_2_dX2
                    + jac2(1, 2)
                          * (jac2(2, 2) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 2) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][1][2]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                  + djac1(3, 1, 3) * l3_1
                  + jac1(1, 1)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 1)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 1)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                    + djac1(3, 2, 3) * l3_1
                    + jac1(1, 2)
                          * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1
                             + jac1(1, 3) * dl1_1_dX1)
                    + jac1(2, 2)
                          * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1
                             + jac1(1, 3) * dl2_1_dX1)
                    + jac1(3, 2)
                          * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                             + jac1(1, 3) * dl3_1_dX1))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][1][3]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][1][3]
      = 2
        * (H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
               * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                  + djac1(3, 1, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 1, 1) + djac1(1, 1, 1) * l1_1 + djac1(2, 1, 1) * l2_1
                    + djac1(3, 1, 1) * l3_1 + Power(jac1(3, 1), 2) * dl3_1_dZ1
                    + Power(jac1(2, 1), 2) * dl2_1_dY1
                    + jac1(2, 1) * jac1(3, 1) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 1), 2) * dl1_1_dX1
                    + jac1(1, 1)
                          * (jac1(2, 1) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 1) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 1, 1) + djac2(1, 1, 1) * l1_2 + djac2(2, 1, 1) * l2_2
                    + djac2(3, 1, 1) * l3_2 + Power(jac2(3, 1), 2) * dl3_2_dZ2
                    + Power(jac2(2, 1), 2) * dl2_2_dY2
                    + jac2(2, 1) * jac2(3, 1) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 1), 2) * dl1_2_dX2
                    + jac2(1, 1)
                          * (jac2(2, 1) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 1) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][1][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                  + djac1(3, 1, 2) * l3_1
                  + jac1(1, 1)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 1)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 1)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                    + djac1(3, 2, 3) * l3_1
                    + jac1(1, 3)
                          * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1
                             + jac1(1, 2) * dl1_1_dX1)
                    + jac1(2, 3)
                          * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1
                             + jac1(1, 2) * dl2_1_dX1)
                    + jac1(3, 3)
                          * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                             + jac1(1, 2) * dl3_1_dX1))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][1][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                  + djac1(3, 1, 3) * l3_1
                  + jac1(1, 1)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 1)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 1)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (djac1(0, 3, 3) + djac1(1, 3, 3) * l1_1 + djac1(2, 3, 3) * l2_1
                    + djac1(3, 3, 3) * l3_1 + Power(jac1(3, 3), 2) * dl3_1_dZ1
                    + Power(jac1(2, 3), 2) * dl2_1_dY1
                    + jac1(2, 3) * jac1(3, 3) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 3), 2) * dl1_1_dX1
                    + jac1(1, 3)
                          * (jac1(2, 3) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 3) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 1)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 1)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 1)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (djac2(0, 3, 3) + djac2(1, 3, 3) * l1_2 + djac2(2, 3, 3) * l2_2
                    + djac2(3, 3, 3) * l3_2 + Power(jac2(3, 3), 2) * dl3_2_dZ2
                    + Power(jac2(2, 3), 2) * dl2_2_dY2
                    + jac2(2, 3) * jac2(3, 3) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 3), 2) * dl1_2_dX2
                    + jac2(1, 3)
                          * (jac2(2, 3) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 3) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 1) + jac1(1, 1) * l1_1 + jac1(2, 1) * l2_1 + jac1(3, 1) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 1) + jac2(1, 1) * l1_2 + jac2(2, 1) * l2_2 + jac2(3, 1) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][2][2]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][2][2]
      = 4 * H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
            * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
               + djac1(3, 1, 2) * l3_1
               + jac1(1, 2)
                     * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
               + jac1(2, 2)
                     * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
               + jac1(3, 2)
                     * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1 + jac1(1, 1) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
              * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                 + djac2(3, 1, 2) * l3_2
                 + jac2(1, 2)
                       * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2 + jac2(1, 1) * dl1_2_dX2)
                 + jac2(2, 2)
                       * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2 + jac2(1, 1) * dl2_2_dX2)
                 + jac2(3, 2)
                       * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2 + jac2(1, 1) * dl3_2_dX2))
        + 2 * Power(jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1, 2)
              * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
        + 2 * Power(jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2, 2)
              * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2);
  dllg[2][2][2]
      = 4 * H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
            * (djac1(0, 2, 2) + djac1(1, 2, 2) * l1_1 + djac1(2, 2, 2) * l2_1
               + djac1(3, 2, 2) * l3_1 + Power(jac1(3, 2), 2) * dl3_1_dZ1
               + Power(jac1(2, 2), 2) * dl2_1_dY1
               + jac1(2, 2) * jac1(3, 2) * (dl2_1_dZ1 + dl3_1_dY1)
               + Power(jac1(1, 2), 2) * dl1_1_dX1
               + jac1(1, 2)
                     * (jac1(2, 2) * (dl1_1_dY1 + dl2_1_dX1)
                        + jac1(3, 2) * (dl1_1_dZ1 + dl3_1_dX1)))
        + 4 * H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
              * (djac2(0, 2, 2) + djac2(1, 2, 2) * l1_2 + djac2(2, 2, 2) * l2_2
                 + djac2(3, 2, 2) * l3_2 + Power(jac2(3, 2), 2) * dl3_2_dZ2
                 + Power(jac2(2, 2), 2) * dl2_2_dY2
                 + jac2(2, 2) * jac2(3, 2) * (dl2_2_dZ2 + dl3_2_dY2)
                 + Power(jac2(1, 2), 2) * dl1_2_dX2
                 + jac2(1, 2)
                       * (jac2(2, 2) * (dl1_2_dY2 + dl2_2_dX2)
                          + jac2(3, 2) * (dl1_2_dZ2 + dl3_2_dX2)))
        + 2 * Power(jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1, 2)
              * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
        + 2 * Power(jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2, 2)
              * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2);
  dllg[3][2][2]
      = 4 * H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
            * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
               + djac1(3, 2, 3) * l3_1
               + jac1(1, 2)
                     * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
               + jac1(2, 2)
                     * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
               + jac1(3, 2)
                     * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1 + jac1(1, 3) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
              * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                 + djac2(3, 2, 3) * l3_2
                 + jac2(1, 2)
                       * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2 + jac2(1, 3) * dl1_2_dX2)
                 + jac2(2, 2)
                       * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2 + jac2(1, 3) * dl2_2_dX2)
                 + jac2(3, 2)
                       * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2 + jac2(1, 3) * dl3_2_dX2))
        + 2 * Power(jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1, 2)
              * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
        + 2 * Power(jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2, 2)
              * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2);

  dllg[0][2][3]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][2][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 1, 2) + djac1(1, 1, 2) * l1_1 + djac1(2, 1, 2) * l2_1
                  + djac1(3, 1, 2) * l3_1
                  + jac1(1, 2)
                        * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
                  + jac1(2, 2)
                        * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
                  + jac1(3, 2)
                        * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                           + jac1(1, 1) * dl3_1_dX1))
           + H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
                    + djac1(3, 1, 3) * l3_1
                    + jac1(1, 3)
                          * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1
                             + jac1(1, 1) * dl1_1_dX1)
                    + jac1(2, 3)
                          * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1
                             + jac1(1, 1) * dl2_1_dX1)
                    + jac1(3, 3)
                          * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1
                             + jac1(1, 1) * dl3_1_dX1))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 1, 2) + djac2(1, 1, 2) * l1_2 + djac2(2, 1, 2) * l2_2
                    + djac2(3, 1, 2) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                    + djac2(3, 1, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2
                             + jac2(1, 1) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2
                             + jac2(1, 1) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2
                             + jac2(1, 1) * dl3_2_dX2))
           + (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
           + (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2));
  dllg[2][2][3]
      = 2
        * (H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
               * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                  + djac1(3, 2, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1
                           + jac1(1, 2) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 2, 2) + djac1(1, 2, 2) * l1_1 + djac1(2, 2, 2) * l2_1
                    + djac1(3, 2, 2) * l3_1 + Power(jac1(3, 2), 2) * dl3_1_dZ1
                    + Power(jac1(2, 2), 2) * dl2_1_dY1
                    + jac1(2, 2) * jac1(3, 2) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 2), 2) * dl1_1_dX1
                    + jac1(1, 2)
                          * (jac1(2, 2) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 2) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2
                             + jac2(1, 2) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2
                             + jac2(1, 2) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2
                             + jac2(1, 2) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 2, 2) + djac2(1, 2, 2) * l1_2 + djac2(2, 2, 2) * l2_2
                    + djac2(3, 2, 2) * l3_2 + Power(jac2(3, 2), 2) * dl3_2_dZ2
                    + Power(jac2(2, 2), 2) * dl2_2_dY2
                    + jac2(2, 2) * jac2(3, 2) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 2), 2) * dl1_2_dX2
                    + jac2(1, 2)
                          * (jac2(2, 2) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 2) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
           + (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2));
  dllg[3][2][3]
      = 2
        * (H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
               * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
                  + djac1(3, 2, 3) * l3_1
                  + jac1(1, 2)
                        * (jac1(3, 3) * dl1_1_dZ1 + jac1(2, 3) * dl1_1_dY1 + jac1(1, 3) * dl1_1_dX1)
                  + jac1(2, 2)
                        * (jac1(3, 3) * dl2_1_dZ1 + jac1(2, 3) * dl2_1_dY1 + jac1(1, 3) * dl2_1_dX1)
                  + jac1(3, 2)
                        * (jac1(3, 3) * dl3_1_dZ1 + jac1(2, 3) * dl3_1_dY1
                           + jac1(1, 3) * dl3_1_dX1))
           + H1 * (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (djac1(0, 3, 3) + djac1(1, 3, 3) * l1_1 + djac1(2, 3, 3) * l2_1
                    + djac1(3, 3, 3) * l3_1 + Power(jac1(3, 3), 2) * dl3_1_dZ1
                    + Power(jac1(2, 3), 2) * dl2_1_dY1
                    + jac1(2, 3) * jac1(3, 3) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 3), 2) * dl1_1_dX1
                    + jac1(1, 3)
                          * (jac1(2, 3) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 3) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                    + djac2(3, 2, 3) * l3_2
                    + jac2(1, 2)
                          * (jac2(3, 3) * dl1_2_dZ2 + jac2(2, 3) * dl1_2_dY2
                             + jac2(1, 3) * dl1_2_dX2)
                    + jac2(2, 2)
                          * (jac2(3, 3) * dl2_2_dZ2 + jac2(2, 3) * dl2_2_dY2
                             + jac2(1, 3) * dl2_2_dX2)
                    + jac2(3, 2)
                          * (jac2(3, 3) * dl3_2_dZ2 + jac2(2, 3) * dl3_2_dY2
                             + jac2(1, 3) * dl3_2_dX2))
           + H2 * (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (djac2(0, 3, 3) + djac2(1, 3, 3) * l1_2 + djac2(2, 3, 3) * l2_2
                    + djac2(3, 3, 3) * l3_2 + Power(jac2(3, 3), 2) * dl3_2_dZ2
                    + Power(jac2(2, 3), 2) * dl2_2_dY2
                    + jac2(2, 3) * jac2(3, 3) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 3), 2) * dl1_2_dX2
                    + jac2(1, 3)
                          * (jac2(2, 3) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 3) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 2) + jac1(1, 2) * l1_1 + jac1(2, 2) * l2_1 + jac1(3, 2) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
           + (jac2(0, 2) + jac2(1, 2) * l1_2 + jac2(2, 2) * l2_2 + jac2(3, 2) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2));

  dllg[0][3][3]
      = 2
        * (H1 * (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
               * (djac1(0, 0, 3) + djac1(1, 0, 3) * l1_1 + djac1(2, 0, 3) * l2_1
                  + djac1(3, 0, 3) * l3_1
                  + jac1(1, 3)
                        * (jac1(3, 0) * dl1_1_dZ1 + jac1(2, 0) * dl1_1_dY1 + jac1(1, 0) * dl1_1_dX1)
                  + jac1(2, 3)
                        * (jac1(3, 0) * dl2_1_dZ1 + jac1(2, 0) * dl2_1_dY1 + jac1(1, 0) * dl2_1_dX1)
                  + jac1(3, 3)
                        * (jac1(3, 0) * dl3_1_dZ1 + jac1(2, 0) * dl3_1_dY1
                           + jac1(1, 0) * dl3_1_dX1))
           + H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (djac1(0, 0, 0) + djac1(1, 0, 0) * l1_1 + djac1(2, 0, 0) * l2_1
                    + djac1(3, 0, 0) * l3_1 + Power(jac1(3, 0), 2) * dl3_1_dZ1
                    + Power(jac1(2, 0), 2) * dl2_1_dY1
                    + jac1(2, 0) * jac1(3, 0) * (dl2_1_dZ1 + dl3_1_dY1)
                    + Power(jac1(1, 0), 2) * dl1_1_dX1
                    + jac1(1, 0)
                          * (jac1(2, 0) * (dl1_1_dY1 + dl2_1_dX1)
                             + jac1(3, 0) * (dl1_1_dZ1 + dl3_1_dX1)))
           + H2 * (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (djac2(0, 0, 3) + djac2(1, 0, 3) * l1_2 + djac2(2, 0, 3) * l2_2
                    + djac2(3, 0, 3) * l3_2
                    + jac2(1, 3)
                          * (jac2(3, 0) * dl1_2_dZ2 + jac2(2, 0) * dl1_2_dY2
                             + jac2(1, 0) * dl1_2_dX2)
                    + jac2(2, 3)
                          * (jac2(3, 0) * dl2_2_dZ2 + jac2(2, 0) * dl2_2_dY2
                             + jac2(1, 0) * dl2_2_dX2)
                    + jac2(3, 3)
                          * (jac2(3, 0) * dl3_2_dZ2 + jac2(2, 0) * dl3_2_dY2
                             + jac2(1, 0) * dl3_2_dX2))
           + H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (djac2(0, 0, 0) + djac2(1, 0, 0) * l1_2 + djac2(2, 0, 0) * l2_2
                    + djac2(3, 0, 0) * l3_2 + Power(jac2(3, 0), 2) * dl3_2_dZ2
                    + Power(jac2(2, 0), 2) * dl2_2_dY2
                    + jac2(2, 0) * jac2(3, 0) * (dl2_2_dZ2 + dl3_2_dY2)
                    + Power(jac2(1, 0), 2) * dl1_2_dX2
                    + jac2(1, 0)
                          * (jac2(2, 0) * (dl1_2_dY2 + dl2_2_dX2)
                             + jac2(3, 0) * (dl1_2_dZ2 + dl3_2_dX2)))
           + (jac1(0, 0) + jac1(1, 0) * l1_1 + jac1(2, 0) * l2_1 + jac1(3, 0) * l3_1)
                 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
                 * (jac1(3, 0) * dH1_dZ1 + jac1(2, 0) * dH1_dY1 + jac1(1, 0) * dH1_dX1)
           + (jac2(0, 0) + jac2(1, 0) * l1_2 + jac2(2, 0) * l2_2 + jac2(3, 0) * l3_2)
                 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
                 * (jac2(3, 0) * dH2_dZ2 + jac2(2, 0) * dH2_dY2 + jac2(1, 0) * dH2_dX2));
  dllg[1][3][3]
      = 4 * H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
            * (djac1(0, 1, 3) + djac1(1, 1, 3) * l1_1 + djac1(2, 1, 3) * l2_1
               + djac1(3, 1, 3) * l3_1
               + jac1(1, 3)
                     * (jac1(3, 1) * dl1_1_dZ1 + jac1(2, 1) * dl1_1_dY1 + jac1(1, 1) * dl1_1_dX1)
               + jac1(2, 3)
                     * (jac1(3, 1) * dl2_1_dZ1 + jac1(2, 1) * dl2_1_dY1 + jac1(1, 1) * dl2_1_dX1)
               + jac1(3, 3)
                     * (jac1(3, 1) * dl3_1_dZ1 + jac1(2, 1) * dl3_1_dY1 + jac1(1, 1) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
              * (djac2(0, 1, 3) + djac2(1, 1, 3) * l1_2 + djac2(2, 1, 3) * l2_2
                 + djac2(3, 1, 3) * l3_2
                 + jac2(1, 3)
                       * (jac2(3, 1) * dl1_2_dZ2 + jac2(2, 1) * dl1_2_dY2 + jac2(1, 1) * dl1_2_dX2)
                 + jac2(2, 3)
                       * (jac2(3, 1) * dl2_2_dZ2 + jac2(2, 1) * dl2_2_dY2 + jac2(1, 1) * dl2_2_dX2)
                 + jac2(3, 3)
                       * (jac2(3, 1) * dl3_2_dZ2 + jac2(2, 1) * dl3_2_dY2 + jac2(1, 1) * dl3_2_dX2))
        + 2 * Power(jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1, 2)
              * (jac1(3, 1) * dH1_dZ1 + jac1(2, 1) * dH1_dY1 + jac1(1, 1) * dH1_dX1)
        + 2 * Power(jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2, 2)
              * (jac2(3, 1) * dH2_dZ2 + jac2(2, 1) * dH2_dY2 + jac2(1, 1) * dH2_dX2);
  dllg[2][3][3]
      = 4 * H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
            * (djac1(0, 2, 3) + djac1(1, 2, 3) * l1_1 + djac1(2, 2, 3) * l2_1
               + djac1(3, 2, 3) * l3_1
               + jac1(1, 3)
                     * (jac1(3, 2) * dl1_1_dZ1 + jac1(2, 2) * dl1_1_dY1 + jac1(1, 2) * dl1_1_dX1)
               + jac1(2, 3)
                     * (jac1(3, 2) * dl2_1_dZ1 + jac1(2, 2) * dl2_1_dY1 + jac1(1, 2) * dl2_1_dX1)
               + jac1(3, 3)
                     * (jac1(3, 2) * dl3_1_dZ1 + jac1(2, 2) * dl3_1_dY1 + jac1(1, 2) * dl3_1_dX1))
        + 4 * H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
              * (djac2(0, 2, 3) + djac2(1, 2, 3) * l1_2 + djac2(2, 2, 3) * l2_2
                 + djac2(3, 2, 3) * l3_2
                 + jac2(1, 3)
                       * (jac2(3, 2) * dl1_2_dZ2 + jac2(2, 2) * dl1_2_dY2 + jac2(1, 2) * dl1_2_dX2)
                 + jac2(2, 3)
                       * (jac2(3, 2) * dl2_2_dZ2 + jac2(2, 2) * dl2_2_dY2 + jac2(1, 2) * dl2_2_dX2)
                 + jac2(3, 3)
                       * (jac2(3, 2) * dl3_2_dZ2 + jac2(2, 2) * dl3_2_dY2 + jac2(1, 2) * dl3_2_dX2))
        + 2 * Power(jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1, 2)
              * (jac1(3, 2) * dH1_dZ1 + jac1(2, 2) * dH1_dY1 + jac1(1, 2) * dH1_dX1)
        + 2 * Power(jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2, 2)
              * (jac2(3, 2) * dH2_dZ2 + jac2(2, 2) * dH2_dY2 + jac2(1, 2) * dH2_dX2);
  dllg[3][3][3]
      = 4 * H1 * (jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1)
            * (djac1(0, 3, 3) + djac1(1, 3, 3) * l1_1 + djac1(2, 3, 3) * l2_1
               + djac1(3, 3, 3) * l3_1 + Power(jac1(3, 3), 2) * dl3_1_dZ1
               + Power(jac1(2, 3), 2) * dl2_1_dY1
               + jac1(2, 3) * jac1(3, 3) * (dl2_1_dZ1 + dl3_1_dY1)
               + Power(jac1(1, 3), 2) * dl1_1_dX1
               + jac1(1, 3)
                     * (jac1(2, 3) * (dl1_1_dY1 + dl2_1_dX1)
                        + jac1(3, 3) * (dl1_1_dZ1 + dl3_1_dX1)))
        + 4 * H2 * (jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2)
              * (djac2(0, 3, 3) + djac2(1, 3, 3) * l1_2 + djac2(2, 3, 3) * l2_2
                 + djac2(3, 3, 3) * l3_2 + Power(jac2(3, 3), 2) * dl3_2_dZ2
                 + Power(jac2(2, 3), 2) * dl2_2_dY2
                 + jac2(2, 3) * jac2(3, 3) * (dl2_2_dZ2 + dl3_2_dY2)
                 + Power(jac2(1, 3), 2) * dl1_2_dX2
                 + jac2(1, 3)
                       * (jac2(2, 3) * (dl1_2_dY2 + dl2_2_dX2)
                          + jac2(3, 3) * (dl1_2_dZ2 + dl3_2_dX2)))
        + 2 * Power(jac1(0, 3) + jac1(1, 3) * l1_1 + jac1(2, 3) * l2_1 + jac1(3, 3) * l3_1, 2)
              * (jac1(3, 3) * dH1_dZ1 + jac1(2, 3) * dH1_dY1 + jac1(1, 3) * dH1_dX1)
        + 2 * Power(jac2(0, 3) + jac2(1, 3) * l1_2 + jac2(2, 3) * l2_2 + jac2(3, 3) * l3_2, 2)
              * (jac2(3, 3) * dH2_dZ2 + jac2(2, 3) * dH2_dY2 + jac2(1, 3) * dH2_dX2);

  dllg[0][1][0] = dllg[0][0][1];
  dllg[1][1][0] = dllg[1][0][1];
  dllg[2][1][0] = dllg[2][0][1];
  dllg[3][1][0] = dllg[3][0][1];

  dllg[0][2][0] = dllg[0][0][2];
  dllg[1][2][0] = dllg[1][0][2];
  dllg[2][2][0] = dllg[2][0][2];
  dllg[3][2][0] = dllg[3][0][2];

  dllg[0][2][1] = dllg[0][1][2];
  dllg[1][2][1] = dllg[1][1][2];
  dllg[2][2][1] = dllg[2][1][2];
  dllg[3][2][1] = dllg[3][1][2];

  dllg[0][3][0] = dllg[0][0][3];
  dllg[1][3][0] = dllg[1][0][3];
  dllg[2][3][0] = dllg[2][0][3];
  dllg[3][3][0] = dllg[3][0][3];

  dllg[0][3][1] = dllg[0][1][3];
  dllg[1][3][1] = dllg[1][1][3];
  dllg[2][3][1] = dllg[2][1][3];
  dllg[3][3][1] = dllg[3][1][3];

  dllg[0][3][2] = dllg[0][2][3];
  dllg[1][3][2] = dllg[1][2][3];
  dllg[2][3][2] = dllg[2][2][3];
  dllg[3][3][2] = dllg[3][2][3];

  return dllg;
}

} // namespace sks_aux

#endif // GRLENSING_SKS_AUX_FUNCTIONS