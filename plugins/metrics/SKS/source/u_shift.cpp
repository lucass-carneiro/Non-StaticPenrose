#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::u_shift(double t, double x, double y, double z)
    -> metric_server::spatial_vector {
  using namespace sks_aux;

  const auto llg = llg_SKS(M1, M2, a1, a2, b, t, x, y, z);
  const auto uug = inverse_symmetric_4x4(llg);

  metric_server::spatial_vector ushift{};

  ushift[0] = -uug[0][1] / uug[0][0];
  ushift[1] = -uug[0][2] / uug[0][0];
  ushift[2] = -uug[0][3] / uug[0][0];

  return ushift;
}