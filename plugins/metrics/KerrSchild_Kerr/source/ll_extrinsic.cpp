#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::ll_extrinsic(double, double x, double y,
                                                                        double z)
    -> metric_server::spatial_matrix {

  using namespace ksk_aux;

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  auto dr_dx = d_r_KS_dx(a, x, y, z);
  auto dr_dy = d_r_KS_dy(a, x, y, z);
  auto dr_dz = d_r_KS_dz(a, x, y, z);

  const auto dH_dx = d_H_KS_dx(M, a, dr_dx, r, z);
  const auto dH_dy = d_H_KS_dy(M, a, dr_dy, r, z);
  const auto dH_dz = d_H_KS_dz(M, a, dr_dz, r, z);

  const auto dl1_dx = d_l1_KS_dx(a, dr_dx, r, x, y);
  const auto dl1_dy = d_l1_KS_dy(a, dr_dy, r, x, y);
  const auto dl1_dz = d_l1_KS_dz(a, dr_dz, r, x, y);

  const auto dl2_dx = d_l2_KS_dx(a, dr_dx, r, x, y);
  const auto dl2_dy = d_l2_KS_dy(a, dr_dy, r, x, y);
  const auto dl2_dz = d_l2_KS_dz(a, dr_dz, r, x, y);

  const auto dl3_dx = d_l3_KS_dx(r, dr_dx, z);
  const auto dl3_dy = d_l3_KS_dy(r, dr_dy, z);
  const auto dl3_dz = d_l3_KS_dz(r, dr_dz, z);

  metric_server::spatial_matrix K{};

  K[0][0] = (2
             * (l1 * dH_dx + H * (Power(l1, 2) * (l3 * dH_dz + l2 * dH_dy + l1 * dH_dx) + dl1_dx)
                + 2 * Power(H, 2) * l1 * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))))
            / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
                   * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));
  K[0][1] = (l1 * dH_dy + l2 * dH_dx
             + H * (dl1_dy + 2 * l1 * l2 * (l3 * dH_dz + l2 * dH_dy + l1 * dH_dx) + dl2_dx)
             + 2 * Power(H, 2)
                   * (Power(l2, 2) * (dl1_dy - dl2_dx)
                      + l1 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx))
                      + l2 * l3 * (dl1_dz - dl3_dx)))
            / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
                   * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));
  K[0][2] = (l1
                 * ((1 + 2 * H * Power(l3, 2)) * dH_dz
                    + 2 * H * l2 * (l3 * dH_dy + H * (-dl2_dz + dl3_dy)))
             + l3 * dH_dx + 2 * Power(H, 2) * l3 * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))
             + H * (dl1_dz + dl3_dx) + 2 * H * Power(l1, 2) * (l3 * dH_dx + H * (-dl1_dz + dl3_dx)))
            / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
                   * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));
  K[1][1] = (2
             * (l2 * dH_dy + H * (dl2_dy + Power(l2, 2) * (l3 * dH_dz + l2 * dH_dy + l1 * dH_dx))
                + 2 * Power(H, 2) * l2 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx))))
            / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
                   * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));
  K[1][2] = (l3 * dH_dy + H * (dl2_dz + dl3_dy)
             + 2 * H * Power(l2, 2) * (l3 * dH_dy + H * (-dl2_dz + dl3_dy))
             + 2 * Power(H, 2) * l3 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx))
             + l2
                   * ((1 + 2 * H * Power(l3, 2)) * dH_dz
                      + 2 * H * l1 * (l3 * dH_dx + H * (-dl1_dz + dl3_dx))))
            / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
                   * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));
  K[2][2]
      = (2
         * (H * Power(l3, 3) * dH_dz + H * dl3_dz + H * Power(l3, 2) * (l2 * dH_dy + l1 * dH_dx)
            + l3 * (dH_dz + 2 * Power(H, 2) * (l2 * (-dl2_dz + dl3_dy) + l1 * (-dl1_dz + dl3_dx)))))
        / Sqrt((1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
               * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))));

  K[1][0] = K[0][1];
  K[2][0] = K[0][2];
  K[2][1] = K[1][2];

  return K;
}