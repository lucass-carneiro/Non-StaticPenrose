#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::l_shift(double t, double x, double y, double z)
    -> metric_server::spatial_vector {

  metric_server::spatial_vector lshift{};

  lshift[0] = llgSKS_01(t, x, y, z);
  lshift[1] = llgSKS_02(t, x, y, z);
  lshift[2] = llgSKS_03(t, x, y, z);

  return lshift;
}