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

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::spatial_christoffel(double, double x,
                                                                               double y, double z)
    -> metric_server::chirstofell_t {

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

  metric_server::chirstofell_t Gamma{};

  Gamma[0][0][0] = (l1
                    * (l1 * dH_dx
                       + 2 * H
                             * (Power(l1, 2) * (l3 * dH_dz + l2 * dH_dy)
                                - l1 * (Power(l2, 2) + Power(l3, 2)) * dH_dx + dl1_dx)
                       + 4 * Power(H, 2) * l1 * (l2 * (dl1_dy - dl2_dx) + l3 * (dl1_dz - dl3_dx))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[0][0][1]
      = (l1
         * (l1 * dH_dy
            + 2 * H
                  * (dl1_dy
                     + l2
                           * (l1 * (l3 * dH_dz + l2 * dH_dy)
                              - (Power(l2, 2) + Power(l3, 2)) * dH_dx))
            + 2 * Power(H, 2)
                  * (l3 * (l1 * (dl2_dz - dl3_dy) + l3 * (dl1_dy - dl2_dx))
                     + 2 * Power(l2, 2) * (dl1_dy - dl2_dx) + l2 * l3 * (dl1_dz - dl3_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[0][0][2]
      = (l1
         * (l1
                * ((1 + 2 * H * Power(l3, 2)) * dH_dz
                   + 2 * H * l2 * (l3 * dH_dy + H * (-dl2_dz + dl3_dy)))
            + 2 * H
                  * ((1 + H * (Power(l2, 2) + 2 * Power(l3, 2))) * dl1_dz
                     - l3 * ((Power(l2, 2) + Power(l3, 2)) * dH_dx + H * l2 * (-dl1_dy + dl2_dx))
                     - H * (Power(l2, 2) + 2 * Power(l3, 2)) * dl3_dx)))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[0][1][1]
      = (l2 * (2 * l1 * dH_dy - l2 * dH_dx)
         + 4 * Power(H, 2) * l2
               * (l1 * l3 * (dl2_dz - dl3_dy) + (Power(l2, 2) + Power(l3, 2)) * (dl1_dy - dl2_dx))
         + 2 * H
               * (l1 * (Power(l2, 2) * l3 * dH_dz + Power(l2, 3) * dH_dy + dl2_dy)
                  - l2 * (-dl1_dy + l2 * (Power(l2, 2) + Power(l3, 2)) * dH_dx + dl2_dx)))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[0][1][2]
      = (l1
             * (l2 * (1 + 2 * H * Power(l3, 2)) * dH_dz + l3 * dH_dy
                + 2 * Power(H, 2) * Power(l3, 2) * (dl2_dz - dl3_dy) + H * (dl2_dz + dl3_dy)
                + 2 * H * Power(l2, 2) * (l3 * dH_dy + H * (-dl2_dz + dl3_dy)))
         + (1 + 2 * H * (Power(l2, 2) + Power(l3, 2)))
               * (-(l2 * l3 * dH_dx) + H * (l3 * (dl1_dy - dl2_dx) + l2 * (dl1_dz - dl3_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[0][2][2] = (2 * l1
                        * (H * Power(l3, 3) * dH_dz + H * dl3_dz + H * l2 * Power(l3, 2) * dH_dy
                           + l3 * (dH_dz + 2 * Power(H, 2) * l2 * (-dl2_dz + dl3_dy)))
                    - l3 * (1 + 2 * H * (Power(l2, 2) + Power(l3, 2)))
                          * (l3 * dH_dx + 2 * H * (-dl1_dz + dl3_dx)))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][0][0] = (l1 * (-(l1 * dH_dy) + 2 * l2 * dH_dx)
                    + 2 * H
                          * (-(Power(l1, 4) * dH_dy) + Power(l1, 2) * l3 * (l2 * dH_dz - l3 * dH_dy)
                             + Power(l1, 3) * l2 * dH_dx + l2 * dl1_dx + l1 * (-dl1_dy + dl2_dx))
                    + 4 * Power(H, 2) * l1
                          * (-((Power(l1, 2) + Power(l3, 2)) * (dl1_dy - dl2_dx))
                             + l2 * l3 * (dl1_dz - dl3_dx)))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][0][1] = (l2
                    * (l2 * dH_dx
                       + 2 * H
                             * (l1
                                    * (-((Power(l1, 2) + Power(l3, 2)) * dH_dy)
                                       + l2 * (l3 * dH_dz + l1 * dH_dx))
                                + dl2_dx)
                       + 2 * Power(H, 2)
                             * (l1 * l3 * (dl2_dz - dl3_dy) + 2 * Power(l1, 2) * (-dl1_dy + dl2_dx)
                                + Power(l3, 2) * (-dl1_dy + dl2_dx) + l2 * l3 * (dl1_dz - dl3_dx))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][0][2]
      = (2 * H * Power(l1, 3) * (-(l3 * dH_dy) + H * (dl2_dz - dl3_dy))
         + l1 * (1 + 2 * H * Power(l3, 2)) * (l2 * dH_dz - l3 * dH_dy + H * (dl2_dz - dl3_dy))
         + l2 * l3 * dH_dx
         + 2 * Power(H, 2) * Power(l3, 2) * (l3 * (-dl1_dy + dl2_dx) + l2 * (dl1_dz - dl3_dx))
         + H * (l3 * (-dl1_dy + dl2_dx) + l2 * (dl1_dz + dl3_dx))
         + 2 * H * Power(l1, 2)
               * (l2 * l3 * dH_dx + H * (l3 * (-dl1_dy + dl2_dx) + l2 * (-dl1_dz + dl3_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][1][1] = (l2
                    * (l2 * dH_dy
                       + 2 * H
                             * (-(l2 * (Power(l1, 2) + Power(l3, 2)) * dH_dy) + dl2_dy
                                + Power(l2, 2) * (l3 * dH_dz + l1 * dH_dx))
                       + 4 * Power(H, 2) * l2 * (l3 * (dl2_dz - dl3_dy) + l1 * (-dl1_dy + dl2_dx))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][1][2] = (l2
                    * (2 * H
                           * ((1 + H * (Power(l1, 2) + 2 * Power(l3, 2))) * dl2_dz
                              - Power(l1, 2) * (l3 * dH_dy + H * dl3_dy)
                              - Power(l3, 2) * (l3 * dH_dy + 2 * H * dl3_dy)
                              + H * l1 * l3 * (-dl1_dy + dl2_dx))
                       + l2
                             * ((1 + 2 * H * Power(l3, 2)) * dH_dz
                                + 2 * H * l1 * (l3 * dH_dx + H * (-dl1_dz + dl3_dx)))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[1][2][2] = (-(l3 * (1 + 2 * H * (Power(l1, 2) + Power(l3, 2)))
                      * (l3 * dH_dy + 2 * H * (-dl2_dz + dl3_dy)))
                    + 2 * l2
                          * (H * Power(l3, 3) * dH_dz + H * dl3_dz + H * l1 * Power(l3, 2) * dH_dx
                             + l3 * (dH_dz + 2 * Power(H, 2) * l1 * (-dl1_dz + dl3_dx))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][0][0]
      = (-2 * H * Power(l1, 4) * dH_dz
         + Power(l1, 2) * (-((1 + 2 * H * Power(l2, 2)) * dH_dz) + 2 * H * l2 * l3 * dH_dy)
         + 2 * H * l3 * dl1_dx + 2 * H * Power(l1, 3) * (l3 * dH_dx + 2 * H * (-dl1_dz + dl3_dx))
         + 2 * l1
               * (l3 * dH_dx + H * (-dl1_dz + dl3_dx)
                  + 2 * Power(H, 2) * l2 * (l3 * (dl1_dy - dl2_dx) + l2 * (-dl1_dz + dl3_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][0][1]
      = (-2 * H * Power(l1, 3) * (l2 * dH_dz + H * (dl2_dz - dl3_dy))
         - l1 * (1 + 2 * H * Power(l2, 2)) * (l2 * dH_dz - l3 * dH_dy + H * (dl2_dz - dl3_dy))
         + l2 * l3 * dH_dx
         + 2 * Power(H, 2) * Power(l2, 2) * (l3 * (dl1_dy - dl2_dx) + l2 * (-dl1_dz + dl3_dx))
         + H * (l3 * (dl1_dy + dl2_dx) + l2 * (-dl1_dz + dl3_dx))
         + 2 * H * Power(l1, 2)
               * (l2 * l3 * dH_dx + H * (l3 * (-dl1_dy + dl2_dx) + l2 * (-dl1_dz + dl3_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][0][2] = (l3
                    * (l3 * dH_dx
                       + 2 * H
                             * (-(Power(l1, 3) * dH_dz) + l1 * l2 * (-(l2 * dH_dz) + l3 * dH_dy)
                                + Power(l1, 2) * l3 * dH_dx + dl3_dx)
                       + 2 * Power(H, 2)
                             * (l1 * l2 * (-dl2_dz + dl3_dy) + 2 * Power(l1, 2) * (-dl1_dz + dl3_dx)
                                + l2 * (l3 * (dl1_dy - dl2_dx) + l2 * (-dl1_dz + dl3_dx)))))
                   / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][1][1]
      = (-2 * H * Power(l2, 4) * dH_dz + 2 * H * l3 * dl2_dy
         + 2 * H * Power(l2, 3) * (l3 * dH_dy + 2 * H * (-dl2_dz + dl3_dy))
         + Power(l2, 2) * (-((1 + 2 * H * Power(l1, 2)) * dH_dz) + 2 * H * l1 * l3 * dH_dx)
         + 2 * l2
               * (l3 * dH_dy + H * (-dl2_dz + dl3_dy)
                  + 2 * Power(H, 2) * l1 * (l1 * (-dl2_dz + dl3_dy) + l3 * (-dl1_dy + dl2_dx))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][1][2]
      = (l3
         * (l3 * dH_dy
            - 2 * H
                  * (Power(l1, 2) * (l2 * dH_dz + H * (dl2_dz - dl3_dy))
                     + Power(l2, 2) * (l2 * dH_dz - l3 * dH_dy + 2 * H * (dl2_dz - dl3_dy)) - dl3_dy
                     + l1
                           * (-(l2 * l3 * dH_dx)
                              + H * (l3 * (dl1_dy - dl2_dx) + l2 * (dl1_dz - dl3_dx))))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  Gamma[2][2][2]
      = (l3
         * (2 * H * dl3_dz + 2 * H * Power(l3, 2) * (l2 * dH_dy + l1 * dH_dx)
            + l3
                  * ((1 - 2 * H * (Power(l1, 2) + Power(l2, 2))) * dH_dz
                     + 4 * Power(H, 2) * (l2 * (-dl2_dz + dl3_dy) + l1 * (-dl1_dz + dl3_dx)))))
        / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));

  Gamma[0][1][0] = Gamma[0][0][1];
  Gamma[0][2][0] = Gamma[0][0][2];
  Gamma[0][2][1] = Gamma[0][1][2];
  Gamma[1][1][0] = Gamma[1][0][1];
  Gamma[1][2][0] = Gamma[1][0][2];
  Gamma[1][2][1] = Gamma[1][1][2];
  Gamma[2][1][0] = Gamma[2][0][1];
  Gamma[2][2][0] = Gamma[2][0][2];
  Gamma[2][2][1] = Gamma[2][1][2];

  return Gamma;
}