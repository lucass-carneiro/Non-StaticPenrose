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

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::grad_ushift(double, double x, double y,
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

  metric_server::spatial_matrix grad_ushift{};

  grad_ushift[0][0] = (-4 * Power(H, 2) * Power(l1, 2) * dl1_dx
                       + 2 * H * (1 + 2 * H * (Power(l2, 2) + Power(l3, 2))) * dl1_dx
                       + 2 * l1 * (dH_dx - 4 * Power(H, 2) * (l2 * dl2_dx + l3 * dl3_dx)))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[0][1] = (2
                       * (-2 * Power(H, 2) * Power(l2, 2) * dl2_dx
                          + H * (1 + 2 * H * (Power(l1, 2) + Power(l3, 2))) * dl2_dx
                          + l2 * (dH_dx - 4 * Power(H, 2) * (l1 * dl1_dx + l3 * dl3_dx))))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[0][2]
      = (2 * l3 * (dH_dx - 4 * Power(H, 2) * (l1 * dl1_dx + l2 * dl2_dx))
         + 2 * H * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) - Power(l3, 2))) * dl3_dx)
        / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[1][0] = (-4 * Power(H, 2) * Power(l1, 2) * dl1_dy
                       + 2 * H * (1 + 2 * H * (Power(l2, 2) + Power(l3, 2))) * dl1_dy
                       + 2 * l1 * (dH_dy - 4 * Power(H, 2) * (l2 * dl2_dy + l3 * dl3_dy)))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[1][1] = (2
                       * (-2 * Power(H, 2) * Power(l2, 2) * dl2_dy
                          + H * (1 + 2 * H * (Power(l1, 2) + Power(l3, 2))) * dl2_dy
                          + l2 * (dH_dy - 4 * Power(H, 2) * (l1 * dl1_dy + l3 * dl3_dy))))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[1][2]
      = (2 * l3 * (dH_dy - 4 * Power(H, 2) * (l1 * dl1_dy + l2 * dl2_dy))
         + 2 * H * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) - Power(l3, 2))) * dl3_dy)
        / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[2][0] = (-4 * Power(H, 2) * Power(l1, 2) * dl1_dz
                       + 2 * H * (1 + 2 * H * (Power(l2, 2) + Power(l3, 2))) * dl1_dz
                       + 2 * l1 * (dH_dz - 4 * Power(H, 2) * (l2 * dl2_dz + l3 * dl3_dz)))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[2][1] = (2
                       * (-2 * Power(H, 2) * Power(l2, 2) * dl2_dz
                          + H * (1 + 2 * H * (Power(l1, 2) + Power(l3, 2))) * dl2_dz
                          + l2 * (dH_dz - 4 * Power(H, 2) * (l1 * dl1_dz + l3 * dl3_dz))))
                      / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);
  grad_ushift[2][2]
      = (2 * l3 * (dH_dz - 4 * Power(H, 2) * (l1 * dl1_dz + l2 * dl2_dz))
         + 2 * H * (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) - Power(l3, 2))) * dl3_dz)
        / Power(1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)), 2);

  return grad_ushift;
}