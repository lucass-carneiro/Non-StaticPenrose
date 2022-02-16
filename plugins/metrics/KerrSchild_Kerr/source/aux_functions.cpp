#include "aux_functions.hpp"

auto ksk_aux::r_KS_part_1(double a, double x, double y, double z) -> double {
  return x * x + y * y + z * z - a * a;
}

auto ksk_aux::r_KS_part_2(double part1, double a, double z) -> double {
  return std::sqrt(4 * a * a * z * z + part1 * part1);
}

auto ksk_aux::r_KS_d_part_1_dxi(double xi) -> double { return 2 * xi; }

auto ksk_aux::r_KS_d_part_2_dx_y(double part1, double part2, double dpart1dx_y) -> double {
  return (part1 + dpart1dx_y) / part2;
}

auto ksk_aux::r_KS_d_part_2_dz(double part1, double part2, double dpart1dz, double a, double z)
    -> double {
  return (8 * a * a * z + 2 * part1 * dpart1dz) / (2 * part2); // NOLINT
}

auto ksk_aux::r_KS(double a, double x, double y, double z) -> double {
  const auto part1 = r_KS_part_1(a, x, y, z);
  const auto part2 = r_KS_part_2(part1, a, z);
  return std::sqrt((part1 + part2) / 2);
}

auto ksk_aux::d_r_KS_dx(double a, double x, double y, double z) -> double {
  const auto part1 = r_KS_part_1(a, x, y, z);
  const auto part2 = r_KS_part_2(part1, a, z);
  const auto dpart1dx = r_KS_d_part_1_dxi(x);
  const auto dpart2dx = r_KS_d_part_2_dx_y(part1, part2, dpart1dx);
  return (dpart1dx + dpart2dx) / (2 * std::sqrt(2 * (part1 + part2)));
}

auto ksk_aux::d_r_KS_dy(double a, double x, double y, double z) -> double {
  auto part1 = r_KS_part_1(a, x, y, z);
  auto part2 = r_KS_part_2(part1, a, z);
  auto dpart1dy = r_KS_d_part_1_dxi(y);
  auto dpart2dy = r_KS_d_part_2_dx_y(part1, part2, dpart1dy);
  return (dpart1dy + dpart2dy) / (2 * std::sqrt(2 * (part1 + part2)));
}

auto ksk_aux::d_r_KS_dz(double a, double x, double y, double z) -> double {
  auto part1 = r_KS_part_1(a, x, y, z);
  auto part2 = r_KS_part_2(part1, a, z);
  auto dpart1dz = r_KS_d_part_1_dxi(z);
  auto dpart2dz = r_KS_d_part_2_dz(part1, part2, dpart1dz, a, z);
  return (dpart1dz + dpart2dz) / (2 * std::sqrt(2 * (part1 + part2)));
}