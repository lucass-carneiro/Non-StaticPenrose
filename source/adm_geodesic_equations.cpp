#include "metric_server.hpp"
#include "time_integration.hpp"

#include <array>
#include <cstddef>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

// TODO: Review geodesic eq and adm metrics

using svector = grlensing::callable_array<3, realtype>;

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

  // Contracted vectors
  const auto cv1 = [&](std::size_t i) -> realtype {
    return grad_lapse(0) * uu_smetric(i, 0) + grad_lapse(1) * uu_smetric(i, 1)
           + grad_lapse(2) * uu_smetric(i, 2);
  };

  const auto cv2 = [&](std::size_t i) -> realtype {
    return grad_ushift(0, i) * V(0) + grad_ushift(1, i) * V(1) + grad_ushift(2, i) * V(2);
  };

  const auto cv3 = [&](std::size_t i) -> realtype {
    return lapse
           * (ul_extrinsic(i, 0) * V(0) + ul_extrinsic(i, 1) * V(1) + ul_extrinsic(i, 2) * V(2));
  };

  const auto cv4 = [&](std::size_t i) -> realtype {
    return lapse
           * (spatial_christoffel(i, 0, 0) * Power(V(0), 2)
              + 2 * spatial_christoffel(i, 0, 1) * V(0) * V(1)
              + spatial_christoffel(i, 1, 1) * Power(V(1), 2)
              + 2 * (spatial_christoffel(i, 0, 2) * V(0) + spatial_christoffel(i, 1, 2) * V(1))
                    * V(2)
              + spatial_christoffel(i, 2, 2) * Power(V(2), 2));
  };

  // Contracted scalars
  const auto cs1 = grad_lapse(0) * V(0) + grad_lapse(1) * V(1) + grad_lapse(2) * V(2);
  const auto cs2 = lapse
                   * (ll_extrinsic(0, 0) * Power(V(0), 2)
                      + V(1) * (2 * ll_extrinsic(0, 1) * V(0) + ll_extrinsic(1, 1) * V(1))
                      + (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)) * V(2)
                      + ll_extrinsic(2, 2) * Power(V(2), 2));

  // NOLINTNEXTLINE
  NV_Ith_S(ydot, 0) = -cv1(0) - cv2(0) + 2 * cv3(0) - cv4(0) + V(0) * (cs1 - cs2);
  NV_Ith_S(ydot, 1) = -cv1(1) - cv2(1) + 2 * cv3(1) - cv4(1) + V(1) * (cs1 - cs2);
  NV_Ith_S(ydot, 2) = -cv1(2) - cv2(2) + 2 * cv3(2) - cv4(2) + V(2) * (cs1 - cs2);
  NV_Ith_S(ydot, 3) = lapse * V(0) - u_shift(0);
  NV_Ith_S(ydot, 4) = lapse * V(1) - u_shift(1);
  NV_Ith_S(ydot, 5) = lapse * V(2) - u_shift(2);
  NV_Ith_S(ydot, 6) = El * (cs2 - cs1);

  return 0;
}
