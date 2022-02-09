#include "time_integration.hpp"

#include <cmath>
#include <sundials/sundials_types.h>

using svector = std::array<realtype, 3>;

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T(1) : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T(1) / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

auto grlensing::detect_collisions(realtype, N_Vector y, realtype *cs, void *user_data) noexcept
    -> int {

  const svector X = {NV_Ith_S(y, 3), NV_Ith_S(y, 4), NV_Ith_S(y, 5)}; // NOLINT
  const auto EN = NV_Ith_S(y, 6);                                     // NOLINT

  // Reinterpret user data
  const auto cast_user_data = static_cast<user_data_t *>(user_data);

  const auto &traj_conf = std::get<1>(*cast_user_data);

  // Background colision
  const realtype current_radius = std::sqrt(Power(X[0], 2) + Power(X[1], 2) + Power(X[2], 2));
  cs[0] = traj_conf.background_radius - current_radius; // NOLINT

  // Swallowing by object
  // NOLINTNEXTLINE
  cs[1] = std::abs(EN - traj_conf.initial_EN) / traj_conf.initial_EN - traj_conf.energy_threshold;

  return 0;
}