#ifndef GRLENSING_KERRSCHILD_KERR_AUX_FUNCTIONS_HPP
#define GRLENSING_KERRSCHILD_KERR_AUX_FUNCTIONS_HPP

#include <cmath>

namespace ksk_aux {

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T{1} : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T{1} / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

constexpr auto Sqrt(auto x) { return std::sqrt(x); }

auto r_KS(double a, double x, double y, double z) -> double;

auto d_r_KS_dx(double a, double x, double y, double z) -> double;
auto d_r_KS_dy(double a, double x, double y, double z) -> double;
auto d_r_KS_dz(double a, double x, double y, double z) -> double;

auto H_KS(double M, double a, double r, double z) -> double;

auto d_H_KS_dx(double M, double a, double dr_dx, double r, double z) -> double;
auto d_H_KS_dy(double M, double a, double dr_dy, double r, double z) -> double;
auto d_H_KS_dz(double M, double a, double dr_dz, double r, double z) -> double;

auto l1_KS(double a, double r, double x, double y) -> double;

auto d_l1_KS_dx(double a, double dr_dx, double r, double x, double y) -> double;
auto d_l1_KS_dy(double a, double dr_dy, double r, double x, double y) -> double;
auto d_l1_KS_dz(double a, double dr_dz, double r, double x, double y) -> double;

auto l2_KS(double a, double r, double x, double y) -> double;

auto d_l2_KS_dx(double a, double dr_dx, double r, double x, double y) -> double;
auto d_l2_KS_dy(double a, double dr_dy, double r, double x, double y) -> double;
auto d_l2_KS_dz(double a, double dr_dz, double r, double x, double y) -> double;

auto l3_KS(double r, double z) -> double;

auto d_l3_KS_dx(double r, double dr_dx, double z) -> double;
auto d_l3_KS_dy(double r, double dr_dy, double z) -> double;
auto d_l3_KS_dz(double r, double dr_dz, double z) -> double;

} // namespace ksk_aux

#endif // GRLENSING_KERRSCHILD_KERR_AUX_FUNCTIONS_HPP