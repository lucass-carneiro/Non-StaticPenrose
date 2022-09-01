#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::uu_smetric(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {

  using std::pow;

  metric_server::spatial_matrix uugamma{};

  const double v0{llgSKS_33(t, x, y, z)};
  const double v1{llgSKS_22(t, x, y, z)};
  const double v2{llgSKS_11(t, x, y, z)};
  const double v3{llgSKS_12(t, x, y, z)};
  const double v4{llgSKS_23(t, x, y, z)};
  const double v5{llgSKS_13(t, x, y, z)};
  const double v6{pow(v3, 2)};
  const double v7{pow(v4, 2)};
  const double v8{pow(v5, 2)};
  const double v9{1 / (-(v0 * v1 * v2) - 2 * v3 * v4 * v5 + v0 * v6 + v2 * v7 + v1 * v8)};
  const double v10{(v2 * v4 - v3 * v5) * v9};
  const double v11{(-(v3 * v4) + v1 * v5) * v9};
  const double v12{(v0 * v3 - v4 * v5) * v9};

  uugamma[0][0] = (-(v0 * v1) + v7) * v9;
  uugamma[0][1] = v12;
  uugamma[0][2] = v11;
  uugamma[1][0] = v12;
  uugamma[1][1] = (-(v0 * v2) + v8) * v9;
  uugamma[1][2] = v10;
  uugamma[2][0] = v11;
  uugamma[2][1] = v10;
  uugamma[2][2] = (-(v1 * v2) + v6) * v9;

  return uugamma;
}