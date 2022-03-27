#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::d_H_KS_dx;
using ksk_aux::d_H_KS_dy;
using ksk_aux::d_H_KS_dz;
using ksk_aux::d_l1_KS_dx;
using ksk_aux::d_l1_KS_dy;
using ksk_aux::d_l1_KS_dz;
using ksk_aux::d_l2_KS_dx;
using ksk_aux::d_l2_KS_dy;
using ksk_aux::d_l2_KS_dz;
using ksk_aux::d_l3_KS_dx;
using ksk_aux::d_l3_KS_dy;
using ksk_aux::d_l3_KS_dz;
using ksk_aux::d_r_KS_dx;
using ksk_aux::H_KS;
using ksk_aux::l1_KS;
using ksk_aux::l2_KS;
using ksk_aux::l3_KS;
using ksk_aux::Power;
using ksk_aux::r_KS;
using ksk_aux::Sqrt;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::ul_extrinsic(double, double x, double y,
                                                                        double z)
    -> metric_server::spatial_matrix {

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  auto dr_dx = d_r_KS_dx(a, x, y, z);

  const auto dH_dx = d_H_KS_dx(M, a, dr_dx, r, z);
  const auto dH_dy = d_H_KS_dy(M, a, dr_dx, r, z);
  const auto dH_dz = d_H_KS_dz(M, a, dr_dx, r, z);

  const auto dl1_dx = d_l1_KS_dx(a, dr_dx, r, x, y);
  const auto dl1_dy = d_l1_KS_dy(a, dr_dx, r, x, y);
  const auto dl1_dz = d_l1_KS_dz(a, dr_dx, r, x, y);

  const auto dl2_dx = d_l2_KS_dx(a, dr_dx, r, x, y);
  const auto dl2_dy = d_l2_KS_dy(a, dr_dx, r, x, y);
  const auto dl2_dz = d_l2_KS_dz(a, dr_dx, r, x, y);

  const auto dl3_dx = d_l3_KS_dx(r, dr_dx, z);
  const auto dl3_dy = d_l3_KS_dy(r, dr_dx, z);
  const auto dl3_dz = d_l3_KS_dz(r, dr_dx, z);

  metric_server::spatial_matrix mixedK{};

  mixedK[0][0]
      = (2
         * (l1 * dH_dx + H * (l1 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)) * dH_dx + dl1_dx)
            + Power(H, 2)
                  * (2 * (Power(l2, 2) + Power(l3, 2)) * dl1_dx
                     + l1 * (l2 * (dl1_dy - 3 * dl2_dx) + l3 * (dl1_dz - 3 * dl3_dx)))
            + 2 * Power(H, 3) * l1 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))
                  * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[0][1]
      = (l1 * dH_dy + l2 * dH_dx
         + H * (dl1_dy + 2 * l2 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)) * dH_dx + dl2_dx)
         + 2 * Power(H, 2)
               * (2 * Power(l2, 2) * dl1_dy - 2 * l1 * l3 * dl3_dy
                  + Power(l1, 2) * (-dl1_dy + dl2_dx) + Power(l3, 2) * (dl1_dy + dl2_dx)
                  + l2 * (-2 * l1 * dl2_dy + l3 * (dl1_dz - dl3_dx)))
         + 4 * Power(H, 3) * l2 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))
               * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx)))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[0][2]
      = (l1 * (dH_dz - 4 * Power(H, 2) * (l2 * dl2_dz + l3 * dl3_dz))
         + 2 * H * Power(l1, 2)
               * (l3 * dH_dx
                  + 2 * Power(H, 2) * l3 * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))
                  + H * (-dl1_dz + dl3_dx))
         + (1 + 2 * H * (Power(l2, 2) + Power(l3, 2)))
               * (l3 * dH_dx
                  + 2 * Power(H, 2) * l3 * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))
                  + H * (dl1_dz + dl3_dx)))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[1][0]
      = (l1 * dH_dy + l2 * dH_dx
         - 4 * Power(H, 3) * l1 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))
               * (l3 * (-dl2_dz + dl3_dy) + l1 * (dl1_dy - dl2_dx))
         + H * (2 * l1 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)) * dH_dy + dl1_dy + dl2_dx)
         + 2 * Power(H, 2)
               * (l1 * (l3 * (dl2_dz - dl3_dy) - 2 * l2 * dl1_dx) + Power(l2, 2) * (dl1_dy - dl2_dx)
                  + 2 * Power(l1, 2) * dl2_dx + Power(l3, 2) * (dl1_dy + dl2_dx)
                  - 2 * l2 * l3 * dl3_dx))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[1][1]
      = (2
         * (l2 * dH_dy + H * (l2 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)) * dH_dy + dl2_dy)
            - 2 * Power(H, 3) * l2 * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2))
                  * (l3 * (-dl2_dz + dl3_dy) + l1 * (dl1_dy - dl2_dx))
            + Power(H, 2)
                  * (2 * (Power(l1, 2) + Power(l3, 2)) * dl2_dy
                     + l2 * (l3 * (dl2_dz - 3 * dl3_dy) + l1 * (-3 * dl1_dy + dl2_dx)))))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[1][2]
      = (l2 * (dH_dz - 4 * Power(H, 2) * (l1 * dl1_dz + l3 * dl3_dz))
         + 2 * H * Power(l2, 2)
               * (l3 * dH_dy + H * (-dl2_dz + dl3_dy)
                  + 2 * Power(H, 2) * l3 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx)))
         + (1 + 2 * H * (Power(l1, 2) + Power(l3, 2)))
               * (l3 * dH_dy + H * (dl2_dz + dl3_dy)
                  + 2 * Power(H, 2) * l3 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx))))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[2][0]
      = (2 * H * Power(l1, 3) * (dH_dz + 2 * Power(H, 2) * l2 * (-dl2_dz + dl3_dy)) + l3 * dH_dx
         + l1
               * ((1 + 2 * H * (Power(l2, 2) + Power(l3, 2)))
                      * (dH_dz + 2 * Power(H, 2) * l2 * (-dl2_dz + dl3_dy))
                  - 4 * Power(H, 2) * l3 * dl1_dx)
         + 4 * Power(H, 3) * Power(l1, 4) * (-dl1_dz + dl3_dx) + H * (dl1_dz + dl3_dx)
         + 4 * Power(H, 2) * Power(l1, 2)
               * (-(H * (Power(l2, 2) + Power(l3, 2)) * (dl1_dz - dl3_dx)) + dl3_dx)
         + 2 * Power(H, 2)
               * (-2 * l2 * l3 * dl2_dx + Power(l3, 2) * (dl1_dz - dl3_dx)
                  + Power(l2, 2) * (dl1_dz + dl3_dx)))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[2][1]
      = (l3 * dH_dy + 4 * Power(H, 3) * Power(l2, 4) * (-dl2_dz + dl3_dy) + H * (dl2_dz + dl3_dy)
         + 4 * Power(H, 2) * Power(l2, 2)
               * (-(H * (Power(l1, 2) + Power(l3, 2)) * (dl2_dz - dl3_dy)) + dl3_dy)
         + 2 * Power(H, 2)
               * (-2 * l1 * l3 * dl1_dy + Power(l3, 2) * (dl2_dz - dl3_dy)
                  + Power(l1, 2) * (dl2_dz + dl3_dy))
         + l2
               * ((1 + 2 * H * (Power(l1, 2) + Power(l3, 2))) * dH_dz
                  - 2 * Power(H, 2)
                        * (2 * l3 * dl2_dy + 2 * H * Power(l1, 3) * (dl1_dz - dl3_dx)
                           + l1 * (1 + 2 * H * Power(l3, 2)) * (dl1_dz - dl3_dx)))
         + Power(l2, 3) * (2 * H * dH_dz + 4 * Power(H, 3) * l1 * (-dl1_dz + dl3_dx)))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));
  mixedK[2][2]
      = (2
         * (H * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2))) * dl3_dz
            + H * Power(l3, 3)
                  * (dH_dz + 2 * Power(H, 2) * (l2 * (-dl2_dz + dl3_dy) + l1 * (-dl1_dz + dl3_dx)))
            + l3
                  * ((1 + H * (Power(l1, 2) + Power(l2, 2))) * dH_dz
                     + Power(H, 2)
                           * (2 * H * Power(l1, 2) * l2 * (-dl2_dz + dl3_dy)
                              + l2
                                    * (-((3 + 2 * H * Power(l2, 2)) * dl2_dz)
                                       + (1 + 2 * H * Power(l2, 2)) * dl3_dy)
                              + 2 * H * Power(l1, 3) * (-dl1_dz + dl3_dx)
                              + l1
                                    * (-((3 + 2 * H * Power(l2, 2)) * dl1_dz)
                                       + (1 + 2 * H * Power(l2, 2)) * dl3_dx)))))
        / (Sqrt(1 + 2 * H * (-1 + Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
           * ((1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))
              * Sqrt(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)))));

  return mixedK;
}