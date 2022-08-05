#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_00(double t, double x, double y, double z) const noexcept -> double {

  using namespace sks_aux;
  using std::pow;
  using std::sqrt;

  const double v0 = Z(2, t, x, y, z);
  const double v1 = Y(2, t, x, y, z);
  const double v2 = X(2, t, x, y, z);
  const double v3 = Z(1, t, x, y, z);
  const double v4 = Y(1, t, x, y, z);
  const double v5 = X(1, t, x, y, z);
  const double v6 = dsy_dt(2, t);
  const double v7 = dsy_dt(1, t);
  const double v8 = dsx_dt(1, t);
  const double v9 = dsx_dt(2, t);
  const double v10 = pow(v6, 2);
  const double v11 = pow(v7, 2);
  const double v12 = pow(v8, 2);
  const double v13 = pow(v9, 2);
  const double v14 = l2_KS(a2, v2, v1, v0);
  const double v15 = l2_KS(a1, v5, v4, v3);
  const double v16 = l1_KS(a1, v5, v4, v3);
  const double v17 = l1_KS(a2, v2, v1, v0);
  const double v18 = H_KS(M2, a2, v2, v1, v0);
  const double v19 = H_KS(M1, a1, v5, v4, v3);
  const double v20 = pow(v14, 2);
  const double v21 = pow(v15, 2);
  const double v22 = pow(v16, 2);
  const double v23 = pow(v17, 2);

  return (-1 + v10 + v11 - v10 * v11 + v12 - v10 * v12 + v13 - v11 * v13 - v12 * v13 + v18
          - v11 * v18 - v12 * v18 + v19 - v10 * v19 - v13 * v19 + v10 * v18 * v20
          - v10 * v11 * v18 * v20 - v10 * v12 * v18 * v20 + v11 * v19 * v21 - v10 * v11 * v19 * v21
          - v11 * v13 * v19 * v21 + v12 * v19 * v22 - v10 * v12 * v19 * v22 - v12 * v13 * v19 * v22
          + v13 * v18 * v23 - v11 * v13 * v18 * v23 - v12 * v13 * v18 * v23 - 2 * v14 * v18 * v6
          + 2 * v11 * v14 * v18 * v6 + 2 * v12 * v14 * v18 * v6 - 2 * v15 * v19 * v7
          + 2 * v10 * v15 * v19 * v7 + 2 * v13 * v15 * v19 * v7 - 2 * v16 * v19 * v8
          + 2 * v10 * v16 * v19 * v8 + 2 * v13 * v16 * v19 * v8 + 2 * v15 * v16 * v19 * v7 * v8
          - 2 * v10 * v15 * v16 * v19 * v7 * v8 - 2 * v13 * v15 * v16 * v19 * v7 * v8
          - 2 * v17 * v18 * v9 + 2 * v11 * v17 * v18 * v9 + 2 * v12 * v17 * v18 * v9
          + 2 * v14 * v17 * v18 * v6 * v9 - 2 * v11 * v14 * v17 * v18 * v6 * v9
          - 2 * v12 * v14 * v17 * v18 * v6 * v9)
         / ((-1 + v11 + v12) * (-1 + v10 + v13));
}