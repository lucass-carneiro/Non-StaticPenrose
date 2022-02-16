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

auto r_KS_part_1(double a, double x, double y, double z) -> double;
auto r_KS_part_2(double part1, double a, double z) -> double;

auto r_KS_d_part_1_dxi(double xi) -> double;
auto r_KS_d_part_2_dx_y(double part1, double part2, double dpart1dx_y) -> double;
auto r_KS_d_part_2_dz(double part1, double part2, double dpart1dz, double a, double z) -> double;

auto r_KS(double a, double x, double y, double z) -> double;

auto d_r_KS_dx(double a, double x, double y, double z) -> double;
auto d_r_KS_dy(double a, double x, double y, double z) -> double;
auto d_r_KS_dz(double a, double x, double y, double z) -> double;

} // namespace ksk_aux

#endif // GRLENSING_KERRSCHILD_KERR_AUX_FUNCTIONS_HPP