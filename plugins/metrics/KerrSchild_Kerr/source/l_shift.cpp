#include "KerrSchild_Kerr.hpp"
#include "aux_functions.hpp"

using grlensing::KerrSchild_Kerr;
using ksk_aux::H_KS;
using ksk_aux::l1_KS;
using ksk_aux::l2_KS;
using ksk_aux::l3_KS;
using ksk_aux::r_KS;

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::l_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {

  auto r = r_KS(a, x, y, z);
  auto H = H_KS(M, a, r, z);
  auto l1 = l1_KS(a, r, x, y);
  auto l2 = l2_KS(a, r, x, y);
  auto l3 = l3_KS(r, z);

  metric_server::spatial_vector lshift{};

  lshift[0] = 2 * H * l1;
  lshift[1] = 2 * H * l2;
  lshift[2] = 2 * H * l3;

  return lshift;
}