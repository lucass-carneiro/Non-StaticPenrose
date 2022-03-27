#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::H_KS;
using ksk_aux::l1_KS;
using ksk_aux::l2_KS;
using ksk_aux::l3_KS;
using ksk_aux::Power;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::ll_smetric(double, double x, double y,
                                                                      double z)
    -> metric_server::spatial_matrix {

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  metric_server::spatial_matrix llsmetric{};

  llsmetric[0][0] = 1 + 2 * H * Power(l1, 2);
  llsmetric[0][1] = 2 * H * l1 * l2;
  llsmetric[0][2] = 2 * H * l1 * l3;
  llsmetric[1][1] = 1 + 2 * H * Power(l2, 2);
  llsmetric[1][2] = 2 * H * l2 * l3;
  llsmetric[2][2] = 1 + 2 * H * Power(l3, 2);

  llsmetric[1][0] = llsmetric[0][1];
  llsmetric[2][0] = llsmetric[0][2];
  llsmetric[2][1] = llsmetric[1][2];

  return llsmetric;
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::uu_smetric(double, double x, double y,
                                                                      double z)
    -> metric_server::spatial_matrix {

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  metric_server::spatial_matrix uusmetric{};

  uusmetric[0][0] = (1 + 2 * H * (Power(l2, 2) + Power(l3, 2)))
                    / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  uusmetric[0][1] = (-2 * H * l1 * l2) / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  uusmetric[0][2] = (-2 * H * l1 * l3) / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  uusmetric[1][1] = (1 + 2 * H * (Power(l1, 2) + Power(l3, 2)))
                    / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  uusmetric[1][2] = (-2 * H * l2 * l3) / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  uusmetric[2][2] = (1 + 2 * H * (Power(l1, 2) + Power(l2, 2)))
                    / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));

  uusmetric[1][0] = uusmetric[0][1];
  uusmetric[2][0] = uusmetric[0][2];
  uusmetric[2][1] = uusmetric[1][2];

  return uusmetric;
}