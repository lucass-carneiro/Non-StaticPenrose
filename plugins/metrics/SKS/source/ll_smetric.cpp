#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::ll_smetric(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {

  metric_server::spatial_matrix llgamma{};

  const double v0 = llgSKS_23(t, x, y, z);
  const double v1 = llgSKS_13(t, x, y, z);
  const double v2 = llgSKS_12(t, x, y, z);

  llgamma[0][0] = llgSKS_11(t, x, y, z);
  llgamma[0][1] = v2;
  llgamma[0][2] = v1;
  llgamma[1][0] = v2;
  llgamma[1][1] = llgSKS_22(t, x, y, z);
  llgamma[1][2] = v0;
  llgamma[2][0] = v1;
  llgamma[2][1] = v0;
  llgamma[2][2] = llgSKS_33(t, x, y, z);

  return llgamma;
}