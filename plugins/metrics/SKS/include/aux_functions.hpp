#ifndef GRLENSING_SKS_AUX_FUNCTIONS
#define GRLENSING_SKS_AUX_FUNCTIONS

#include "../../../../include/tensor_types.hpp"

#include <cmath>
#include <cstddef>

namespace sks_aux {

enum class bh_tag : int { A = 1, B = 2 };

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T{1} : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T{1} / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

constexpr auto Sqrt(auto x) { return std::sqrt(x); }

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

  grlensing::callable_matrix<4, 4, double> jac{};

  jac[0][0] = (4 * sqrt(4 - Power(b, 2) * Power(Omega, 2))
               + Power(-1, static_cast<int>(tag)) * b * x * Power(Omega, 2)
                     * (2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(t * Omega)
               - 3 * Power(-1, static_cast<int>(tag)) * b * x * Power(Omega, 2)
                     * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(3 * t * Omega)
               + b * Power(Omega, 2)
                     * (-2 * b * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(4 * t * Omega)
                        + 2 * Power(-1, static_cast<int>(tag)) * y
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
        / 8;

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
                  / 4;

  djac[0][2][0] = (Power(-1, static_cast<int>(tag)) * b * Power(Omega, 2)
                   * (2 * (-1 + sqrt(4 - Power(b, 2) * Power(Omega, 2)))
                      + 3 * (-2 + sqrt(4 - Power(b, 2) * Power(Omega, 2))) * cos(2 * t * Omega))
                   * sin(t * Omega))
                  / 4;

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

} // namespace sks_aux

#endif // GRLENSING_SKS_AUX_FUNCTIONS