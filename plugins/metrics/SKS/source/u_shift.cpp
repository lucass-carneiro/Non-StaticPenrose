#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::u_shift(double t, double x, double y, double z)
    -> metric_server::spatial_vector {

  using std::pow;

  metric_server::spatial_vector ushift{};

  const double v0 = llgSKS_33(t, x, y, z);
  const double v1 = llgSKS_22(t, x, y, z);
  const double v2 = llgSKS_11(t, x, y, z);
  const double v3 = llgSKS_12(t, x, y, z);
  const double v4 = llgSKS_23(t, x, y, z);
  const double v5 = llgSKS_13(t, x, y, z);
  const double v6 = llgSKS_01(t, x, y, z);
  const double v7 = llgSKS_02(t, x, y, z);
  const double v8 = llgSKS_03(t, x, y, z);
  const double v9 = pow(v3, 2);
  const double v10 = pow(v4, 2);
  const double v11 = pow(v5, 2);
  const double v12 = 1 / (v1 * v11 - v0 * v1 * v2 + v10 * v2 - 2 * v3 * v4 * v5 + v0 * v9);

  ushift[0]
      = v12
        * (-(v0 * v1 * v6) + v10 * v6 + v0 * v3 * v7 - v4 * v5 * v7 - v3 * v4 * v8 + v1 * v5 * v8);
  ushift[1]
      = v12 * (v0 * v3 * v6 - v4 * v5 * v6 + v11 * v7 - v0 * v2 * v7 + v2 * v4 * v8 - v3 * v5 * v8);
  ushift[2]
      = v12
        * (-(v3 * v4 * v6) + v1 * v5 * v6 + v2 * v4 * v7 - v3 * v5 * v7 - v1 * v2 * v8 + v8 * v9);

  return ushift;
}