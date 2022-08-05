#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::lapse(double t, double x, double y, double z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = llgSKS_33(t, x, y, z);
  const double v1 = llgSKS_22(t, x, y, z);
  const double v2 = llgSKS_11(t, x, y, z);
  const double v3 = llgSKS_00(t, x, y, z);
  const double v4 = llgSKS_01(t, x, y, z);
  const double v5 = llgSKS_12(t, x, y, z);
  const double v6 = llgSKS_02(t, x, y, z);
  const double v7 = llgSKS_23(t, x, y, z);
  const double v8 = llgSKS_13(t, x, y, z);
  const double v9 = llgSKS_03(t, x, y, z);
  const double v10 = pow(v4, 2);
  const double v11 = pow(v5, 2);
  const double v12 = pow(v6, 2);
  const double v13 = pow(v7, 2);
  const double v14 = pow(v8, 2);
  const double v15 = pow(v9, 2);

  return sqrt((-(v0 * v1 * v10) + v10 * v13 + v12 * v14 + v11 * v15 - v0 * v12 * v2 - v1 * v15 * v2
               - v0 * v11 * v3 - v1 * v14 * v3 + v0 * v1 * v2 * v3 - v13 * v2 * v3
               + 2 * v0 * v4 * v5 * v6 + 2 * v3 * v5 * v7 * v8 - 2 * v4 * v6 * v7 * v8
               - 2 * v4 * v5 * v7 * v9 + 2 * v2 * v6 * v7 * v9 + 2 * v1 * v4 * v8 * v9
               - 2 * v5 * v6 * v8 * v9)
              / (v0 * v11 + v1 * v14 - v0 * v1 * v2 + v13 * v2 - 2 * v5 * v7 * v8));
}