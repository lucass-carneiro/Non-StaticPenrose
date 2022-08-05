#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_03(double t, double x, double y, double z) const noexcept -> double {

  using namespace sks_aux;
  using std::pow;
  using std::sqrt;

  const double v0 = Z(1, t, x, y, z);
  const double v1 = Y(1, t, x, y, z);
  const double v2 = X(1, t, x, y, z);
  const double v3 = Z(2, t, x, y, z);
  const double v4 = Y(2, t, x, y, z);
  const double v5 = X(2, t, x, y, z);
  const double v6 = dsy_dt(2, t);
  const double v7 = dsx_dt(2, t);
  const double v8 = dsy_dt(1, t);
  const double v9 = dsx_dt(1, t);
  const double v10 = l3_KS(a1, v2, v1, v0);
  const double v11 = l3_KS(a2, v5, v4, v3);
  const double v12 = H_KS(M1, a1, v2, v1, v0);
  const double v13 = H_KS(M2, a2, v5, v4, v3);
  const double v14 = 1 - pow(v6, 2) - pow(v7, 2);
  const double v15 = 1 - pow(v8, 2) - pow(v9, 2);
  const double v16 = sqrt(v14);
  const double v17 = sqrt(v15);

  return (v10 * v12 * v16 + v11 * v13 * v17 - v10 * v12 * v16 * v9 * l1_KS(a1, v2, v1, v0)
          - v11 * v13 * v17 * v7 * l1_KS(a2, v5, v4, v3)
          - v10 * v12 * v16 * v8 * l2_KS(a1, v2, v1, v0)
          - v11 * v13 * v17 * v6 * l2_KS(a2, v5, v4, v3))
         / (sqrt(v14) * sqrt(v15));
}