#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using ksk_aux::H_KS;
using ksk_aux::l1_KS;
using ksk_aux::l2_KS;
using ksk_aux::l3_KS;
using ksk_aux::Power;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::u_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  metric_server::spatial_vector ushift{};

  const auto constant_factor = (2 * H) / (1 + 2 * H * (Power(l1, 2) + Power(l2, 2) + Power(l3, 2)));
  ushift[0] = constant_factor * l1;
  ushift[1] = constant_factor * l2;
  ushift[2] = constant_factor * l3;

  return ushift;
}