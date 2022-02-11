#include "metric_server.hpp"
#include "time_integration.hpp"

#include <array>
#include <cstddef>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

using svector = std::array<realtype, 3>;

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T(1) : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T(1) / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

auto grlensing::adm_geodesic_system(realtype t, N_Vector y, N_Vector ydot, void *user_data) noexcept
    -> int {

  const svector V = {NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2)}; // NOLINT
  const svector X = {NV_Ith_S(y, 3), NV_Ith_S(y, 4), NV_Ith_S(y, 5)}; // NOLINT
  const auto El = NV_Ith_S(y, 6);                                     // NOLINT

  // Reinterpret user data
  const auto cast_user_data = static_cast<user_data_t *>(user_data);

  const auto &metric = std::get<0>(*cast_user_data);
  const auto lapse = metric->lapse(t, X[0], X[1], X[2]);
  const auto u_shift = metric->u_shift(t, X[0], X[1], X[2]);
  const auto uu_smetric = metric->uu_smetric(t, X[0], X[1], X[2]);
  const auto ll_extrinsic = metric->ll_extrinsic(t, X[0], X[1], X[2]);
  const auto ul_extrinsic = metric->ul_extrinsic(t, X[0], X[1], X[2]);

  const auto grad_lapse = metric->grad_lapse(t, X[0], X[1], X[2]);
  const auto grad_ushift = metric->grad_ushift(t, X[0], X[1], X[2]);
  const auto spatial_christoffel = metric->spatial_christoffel(t, X[0], X[1], X[2]);

  // NOLINTNEXTLINE
  NV_Ith_S(ydot, 0)
      = -(grad_lapse[1] * uu_smetric[1][0]) - grad_lapse[2] * uu_smetric[2][0]
        - V[0]
              * (grad_ushift[0][0] - 2 * lapse * ul_extrinsic[0][0]
                 + lapse * spatial_christoffel[0][0][0] * V[0])
        + grad_lapse[0] * (-uu_smetric[0][0] + Power(V[0], 2)) - grad_ushift[1][0] * V[1]
        - grad_ushift[2][0] * V[2] + V[0] * (grad_lapse[1] * V[1] + grad_lapse[2] * V[2])
        + lapse
              * (-(ll_extrinsic[0][0] * Power(V[0], 3))
                 - V[1]
                       * (-2 * ul_extrinsic[0][1]
                          + 2 * V[0] * (spatial_christoffel[0][1][0] + ll_extrinsic[1][0] * V[0])
                          + (spatial_christoffel[0][1][1] + ll_extrinsic[1][1] * V[0]) * V[1])
                 - 2
                       * (-ul_extrinsic[0][2] + spatial_christoffel[0][2][1] * V[1]
                          + V[0]
                                * (spatial_christoffel[0][2][0] + ll_extrinsic[2][0] * V[0]
                                   + ll_extrinsic[2][1] * V[1]))
                       * V[2]
                 - (spatial_christoffel[0][2][2] + ll_extrinsic[2][2] * V[0]) * Power(V[2], 2));

  // NOLINTNEXTLINE
  NV_Ith_S(ydot, 1)
      = -(grad_lapse[1] * uu_smetric[1][1]) - grad_lapse[2] * uu_smetric[2][1]
        - V[0]
              * (grad_ushift[0][1] - 2 * lapse * ul_extrinsic[1][0]
                 + lapse * spatial_christoffel[1][0][0] * V[0])
        - grad_ushift[1][1] * V[1]
        + 2 * lapse * (ul_extrinsic[1][1] - spatial_christoffel[1][1][0] * V[0]) * V[1]
        + grad_lapse[0] * (-uu_smetric[1][0] + V[0] * V[1])
        - V[1]
              * (lapse * ll_extrinsic[0][0] * Power(V[0], 2) - grad_lapse[1] * V[1]
                 + lapse * V[1]
                       * (spatial_christoffel[1][1][1] + 2 * ll_extrinsic[1][0] * V[0]
                          + ll_extrinsic[1][1] * V[1]))
        - grad_ushift[2][1] * V[2]
        + (grad_lapse[2] * V[1]
           - 2 * lapse
                 * (-ul_extrinsic[1][2] + spatial_christoffel[1][2][0] * V[0]
                    + V[1]
                          * (spatial_christoffel[1][2][1] + ll_extrinsic[2][0] * V[0]
                             + ll_extrinsic[2][1] * V[1])))
              * V[2]
        - lapse * (spatial_christoffel[1][2][2] + ll_extrinsic[2][2] * V[1]) * Power(V[2], 2);

  // NOLINTNEXTLINE
  NV_Ith_S(ydot, 2)
      = -(grad_lapse[1] * uu_smetric[2][1]) - grad_lapse[2] * uu_smetric[2][2]
        - grad_ushift[0][2] * V[0] - grad_ushift[1][2] * V[1] - grad_ushift[2][2] * V[2]
        + grad_lapse[0] * (-uu_smetric[2][0] + V[0] * V[2])
        + lapse
              * (2 * ul_extrinsic[2][0] * V[0] - spatial_christoffel[2][0][0] * Power(V[0], 2)
                 - V[1]
                       * (-2 * ul_extrinsic[2][1] + 2 * spatial_christoffel[2][1][0] * V[0]
                          + spatial_christoffel[2][1][1] * V[1])
                 + 2 * (ul_extrinsic[2][2] - spatial_christoffel[2][2][0] * V[0]) * V[2])
        + V[2]
              * (grad_lapse[1] * V[1] + grad_lapse[2] * V[2]
                 - lapse
                       * (ll_extrinsic[0][0] * Power(V[0], 2)
                          + V[1]
                                * (2 * spatial_christoffel[2][2][1] + 2 * ll_extrinsic[1][0] * V[0]
                                   + ll_extrinsic[1][1] * V[1])
                          + (spatial_christoffel[2][2][2] + 2 * ll_extrinsic[2][0] * V[0]
                             + 2 * ll_extrinsic[2][1] * V[1])
                                * V[2]
                          + ll_extrinsic[2][2] * Power(V[2], 2)));

  NV_Ith_S(ydot, 3) = -u_shift[0] + lapse * V[0]; // NOLINT

  NV_Ith_S(ydot, 4) = -u_shift[1] + lapse * V[1]; // NOLINT

  NV_Ith_S(ydot, 5) = -u_shift[2] + lapse * V[2]; // NOLINT

  // NOLINTNEXTLINE
  NV_Ith_S(ydot, 6)
      = El
        * (-(grad_lapse[0] * V[0]) - grad_lapse[1] * V[1] - grad_lapse[2] * V[2]
           + lapse
                 * (ll_extrinsic[0][0] * Power(V[0], 2)
                    + V[1] * (2 * ll_extrinsic[1][0] * V[0] + ll_extrinsic[1][1] * V[1])
                    + (2 * ll_extrinsic[2][0] * V[0] + 2 * ll_extrinsic[2][1] * V[1]) * V[2]
                    + ll_extrinsic[2][2] * Power(V[2], 2)));

  return 0;
}
