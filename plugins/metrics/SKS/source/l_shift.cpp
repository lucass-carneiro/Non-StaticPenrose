#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::l_shift(double t, double x, double y, double z)
    -> metric_server::spatial_vector {
  using namespace sks_aux;

  const auto llg = llg_SKS(M1, M2, a1, a2, b, t, x, y, z);

  metric_server::spatial_vector lshift{};

  lshift[0] = llg[0][1];
  lshift[1] = llg[0][2];
  lshift[2] = llg[0][3];

  return lshift;
}