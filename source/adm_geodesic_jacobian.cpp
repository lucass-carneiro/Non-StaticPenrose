#include "time_integration.hpp"

using svector = grlensing::callable_array<3, realtype>;

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T(1) : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T(1) / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

auto grlensing::adm_geodesic_jacobian(realtype t, N_Vector y, N_Vector, SUNMatrix J,
                                      void *user_data, N_Vector, N_Vector, N_Vector) noexcept
    -> int {

  const svector V = {NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2)}; // NOLINT
  const svector X = {NV_Ith_S(y, 3), NV_Ith_S(y, 4), NV_Ith_S(y, 5)}; // NOLINT
  const auto EN = NV_Ith_S(y, 6);                                     // NOLINT

  // Reinterpret user data
  const auto cast_user_data = static_cast<user_data_t *>(user_data);

  const auto &metric = std::get<0>(*cast_user_data);
  const auto lapse = metric->lapse(t, X[0], X[1], X[2]);
  const auto ll_extrinsic = metric->ll_extrinsic(t, X[0], X[1], X[2]);
  const auto ul_extrinsic = metric->ul_extrinsic(t, X[0], X[1], X[2]);

  const auto grad_lapse = metric->grad_lapse(t, X[0], X[1], X[2]);
  const auto grad_ushift = metric->grad_ushift(t, X[0], X[1], X[2]);
  const auto spatial_christoffel = metric->spatial_christoffel(t, X[0], X[1], X[2]);

  // Zero fill the Jacobian
  SUNMatZero(J);

  SM_ELEMENT_D(J, 0, 0)
      = -grad_ushift(0, 0) + 2 * lapse * ul_extrinsic(0, 0) + grad_lapse(0) * V(0)
        + grad_lapse(1) * V(1) + grad_lapse(2) * V(2)
        - lapse
              * (2 * spatial_christoffel(0, 0, 0) * V(0) + 2 * spatial_christoffel(0, 0, 1) * V(1)
                 + 2 * spatial_christoffel(0, 0, 2) * V(2))
        - lapse
              * (ll_extrinsic(0, 0) * Power(V(0), 2)
                 + V(1) * (2 * ll_extrinsic(0, 1) * V(0) + ll_extrinsic(1, 1) * V(1))
                 + (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)) * V(2)
                 + ll_extrinsic(2, 2) * Power(V(2), 2))
        + V(0)
              * (grad_lapse(0)
                 - lapse
                       * (2 * ll_extrinsic(0, 0) * V(0) + 2 * ll_extrinsic(0, 1) * V(1)
                          + 2 * ll_extrinsic(0, 2) * V(2)));

  SM_ELEMENT_D(J, 0, 1)
      = -grad_ushift(1, 0) + 2 * lapse * ul_extrinsic(0, 1)
        - lapse
              * (2 * spatial_christoffel(0, 0, 1) * V(0) + 2 * spatial_christoffel(0, 1, 1) * V(1)
                 + 2 * spatial_christoffel(0, 1, 2) * V(2))
        + V(0)
              * (grad_lapse(1)
                 - lapse
                       * (2 * ll_extrinsic(0, 1) * V(0) + 2 * ll_extrinsic(1, 1) * V(1)
                          + 2 * ll_extrinsic(1, 2) * V(2)));

  SM_ELEMENT_D(J, 0, 2)
      = -grad_ushift(2, 0) + 2 * lapse * ul_extrinsic(0, 2)
        - lapse
              * (2 * (spatial_christoffel(0, 0, 2) * V(0) + spatial_christoffel(0, 1, 2) * V(1))
                 + 2 * spatial_christoffel(0, 2, 2) * V(2))
        + V(0)
              * (grad_lapse(2)
                 - lapse
                       * (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)
                          + 2 * ll_extrinsic(2, 2) * V(2)));

  SM_ELEMENT_D(J, 0, 3) = 0;

  SM_ELEMENT_D(J, 0, 4) = 0;

  SM_ELEMENT_D(J, 0, 5) = 0;

  SM_ELEMENT_D(J, 0, 6) = 0;

  SM_ELEMENT_D(J, 1, 0)
      = -grad_ushift(0, 1) + 2 * lapse * ul_extrinsic(1, 0)
        - lapse
              * (2 * spatial_christoffel(1, 0, 0) * V(0) + 2 * spatial_christoffel(1, 0, 1) * V(1)
                 + 2 * spatial_christoffel(1, 0, 2) * V(2))
        + V(1)
              * (grad_lapse(0)
                 - lapse
                       * (2 * ll_extrinsic(0, 0) * V(0) + 2 * ll_extrinsic(0, 1) * V(1)
                          + 2 * ll_extrinsic(0, 2) * V(2)));

  SM_ELEMENT_D(J, 1, 1)
      = -grad_ushift(1, 1) + 2 * lapse * ul_extrinsic(1, 1) + grad_lapse(0) * V(0)
        + grad_lapse(1) * V(1) + grad_lapse(2) * V(2)
        - lapse
              * (2 * spatial_christoffel(1, 0, 1) * V(0) + 2 * spatial_christoffel(1, 1, 1) * V(1)
                 + 2 * spatial_christoffel(1, 1, 2) * V(2))
        - lapse
              * (ll_extrinsic(0, 0) * Power(V(0), 2)
                 + V(1) * (2 * ll_extrinsic(0, 1) * V(0) + ll_extrinsic(1, 1) * V(1))
                 + (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)) * V(2)
                 + ll_extrinsic(2, 2) * Power(V(2), 2))
        + V(1)
              * (grad_lapse(1)
                 - lapse
                       * (2 * ll_extrinsic(0, 1) * V(0) + 2 * ll_extrinsic(1, 1) * V(1)
                          + 2 * ll_extrinsic(1, 2) * V(2)));

  SM_ELEMENT_D(J, 1, 2)
      = -grad_ushift(2, 1) + 2 * lapse * ul_extrinsic(1, 2)
        - lapse
              * (2 * (spatial_christoffel(1, 0, 2) * V(0) + spatial_christoffel(1, 1, 2) * V(1))
                 + 2 * spatial_christoffel(1, 2, 2) * V(2))
        + V(1)
              * (grad_lapse(2)
                 - lapse
                       * (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)
                          + 2 * ll_extrinsic(2, 2) * V(2)));

  SM_ELEMENT_D(J, 1, 3) = 0;

  SM_ELEMENT_D(J, 1, 4) = 0;

  SM_ELEMENT_D(J, 1, 5) = 0;

  SM_ELEMENT_D(J, 1, 6) = 0;

  SM_ELEMENT_D(J, 2, 0)
      = -grad_ushift(0, 2) + 2 * lapse * ul_extrinsic(2, 0)
        - lapse
              * (2 * spatial_christoffel(2, 0, 0) * V(0) + 2 * spatial_christoffel(2, 0, 1) * V(1)
                 + 2 * spatial_christoffel(2, 0, 2) * V(2))
        + V(2)
              * (grad_lapse(0)
                 - lapse
                       * (2 * ll_extrinsic(0, 0) * V(0) + 2 * ll_extrinsic(0, 1) * V(1)
                          + 2 * ll_extrinsic(0, 2) * V(2)));

  SM_ELEMENT_D(J, 2, 1)
      = -grad_ushift(1, 2) + 2 * lapse * ul_extrinsic(2, 1)
        - lapse
              * (2 * spatial_christoffel(2, 0, 1) * V(0) + 2 * spatial_christoffel(2, 1, 1) * V(1)
                 + 2 * spatial_christoffel(2, 1, 2) * V(2))
        + V(2)
              * (grad_lapse(1)
                 - lapse
                       * (2 * ll_extrinsic(0, 1) * V(0) + 2 * ll_extrinsic(1, 1) * V(1)
                          + 2 * ll_extrinsic(1, 2) * V(2)));

  SM_ELEMENT_D(J, 2, 2)
      = -grad_ushift(2, 2) + 2 * lapse * ul_extrinsic(2, 2) + grad_lapse(0) * V(0)
        + grad_lapse(1) * V(1) + grad_lapse(2) * V(2)
        - lapse
              * (2 * (spatial_christoffel(2, 0, 2) * V(0) + spatial_christoffel(2, 1, 2) * V(1))
                 + 2 * spatial_christoffel(2, 2, 2) * V(2))
        - lapse
              * (ll_extrinsic(0, 0) * Power(V(0), 2)
                 + V(1) * (2 * ll_extrinsic(0, 1) * V(0) + ll_extrinsic(1, 1) * V(1))
                 + (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)) * V(2)
                 + ll_extrinsic(2, 2) * Power(V(2), 2))
        + V(2)
              * (grad_lapse(2)
                 - lapse
                       * (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)
                          + 2 * ll_extrinsic(2, 2) * V(2)));

  SM_ELEMENT_D(J, 2, 3) = 0;

  SM_ELEMENT_D(J, 2, 4) = 0;

  SM_ELEMENT_D(J, 2, 5) = 0;

  SM_ELEMENT_D(J, 2, 6) = 0;

  SM_ELEMENT_D(J, 3, 0) = lapse;

  SM_ELEMENT_D(J, 3, 1) = 0;

  SM_ELEMENT_D(J, 3, 2) = 0;

  SM_ELEMENT_D(J, 3, 3) = 0;

  SM_ELEMENT_D(J, 3, 4) = 0;

  SM_ELEMENT_D(J, 3, 5) = 0;

  SM_ELEMENT_D(J, 3, 6) = 0;

  SM_ELEMENT_D(J, 4, 0) = 0;

  SM_ELEMENT_D(J, 4, 1) = lapse;

  SM_ELEMENT_D(J, 4, 2) = 0;

  SM_ELEMENT_D(J, 4, 3) = 0;

  SM_ELEMENT_D(J, 4, 4) = 0;

  SM_ELEMENT_D(J, 4, 5) = 0;

  SM_ELEMENT_D(J, 4, 6) = 0;

  SM_ELEMENT_D(J, 5, 0) = 0;

  SM_ELEMENT_D(J, 5, 1) = 0;

  SM_ELEMENT_D(J, 5, 2) = lapse;

  SM_ELEMENT_D(J, 5, 3) = 0;

  SM_ELEMENT_D(J, 5, 4) = 0;

  SM_ELEMENT_D(J, 5, 5) = 0;

  SM_ELEMENT_D(J, 5, 6) = 0;

  SM_ELEMENT_D(J, 6, 0) = EN
                          * (-grad_lapse(0)
                             + lapse
                                   * (2 * ll_extrinsic(0, 0) * V(0) + 2 * ll_extrinsic(0, 1) * V(1)
                                      + 2 * ll_extrinsic(0, 2) * V(2)));

  SM_ELEMENT_D(J, 6, 1) = EN
                          * (-grad_lapse(1)
                             + lapse
                                   * (2 * ll_extrinsic(0, 1) * V(0) + 2 * ll_extrinsic(1, 1) * V(1)
                                      + 2 * ll_extrinsic(1, 2) * V(2)));

  SM_ELEMENT_D(J, 6, 2) = EN
                          * (-grad_lapse(2)
                             + lapse
                                   * (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)
                                      + 2 * ll_extrinsic(2, 2) * V(2)));

  SM_ELEMENT_D(J, 6, 3) = 0;

  SM_ELEMENT_D(J, 6, 4) = 0;

  SM_ELEMENT_D(J, 6, 5) = 0;

  SM_ELEMENT_D(J, 6, 6)
      = -(grad_lapse(0) * V(0)) - grad_lapse(1) * V(1) - grad_lapse(2) * V(2)
        + lapse
              * (ll_extrinsic(0, 0) * Power(V(0), 2)
                 + V(1) * (2 * ll_extrinsic(0, 1) * V(0) + ll_extrinsic(1, 1) * V(1))
                 + (2 * ll_extrinsic(0, 2) * V(0) + 2 * ll_extrinsic(1, 2) * V(1)) * V(2)
                 + ll_extrinsic(2, 2) * Power(V(2), 2));

  return 0;
}