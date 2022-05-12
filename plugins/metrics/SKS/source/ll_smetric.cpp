#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::ll_smetric(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {
  using namespace sks_aux;

  const auto llg = llg_SKS(M1, M2, a1, a2, b, t, x, y, z);

  metric_server::spatial_matrix llgamma{};

  llgamma[0][0] = llg[1][1];
  llgamma[0][1] = llg[1][2];
  llgamma[0][2] = llg[1][3];
  llgamma[1][1] = llg[2][2];
  llgamma[1][2] = llg[2][3];
  llgamma[2][2] = llg[3][3];

  llgamma[1][0] = llgamma[0][1];
  llgamma[2][0] = llgamma[0][2];
  llgamma[2][1] = llgamma[1][2];

  return llgamma;
}