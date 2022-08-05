#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_13(double t, double x, double y, double z) const noexcept -> double {

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
  const double v7 = dsx_dt(2, t);
  const double v8 = dsy_dt(1, t);
  const double v9 = dsx_dt(1, t);
  const double v10 = pow(v6, 2);
  const double v11 = pow(v7, 2);
  const double v12 = pow(v8, 2);
  const double v13 = pow(v9, 2);
  const double v14 = pow(v9, 3);
  const double v15 = pow(v7, 3);
  const double v16 = -v10;
  const double v17 = -v11;
  const double v18 = -v12;
  const double v19 = -v13;
  const double v20 = l3_KS(a2, v2, v1, v0);
  const double v21 = l1_KS(a2, v2, v1, v0);
  const double v22 = l3_KS(a1, v5, v4, v3);
  const double v23 = l1_KS(a1, v5, v4, v3);
  const double v24 = l2_KS(a1, v5, v4, v3);
  const double v25 = l2_KS(a2, v2, v1, v0);
  const double v26 = H_KS(M2, a2, v2, v1, v0);
  const double v27 = H_KS(M1, a1, v5, v4, v3);
  const double v28 = 1 + v16 + v17;
  const double v29 = 1 + v18 + v19;
  const double v30 = sqrt(v28);
  const double v31 = sqrt(v29);

  return (v15 * sqrt(1 - v13 + v18) * v19 * v20 * v26 + v15 * v18 * sqrt(1 - v12 + v19) * v20 * v26
          + v14 * sqrt(1 - v11 + v16) * v17 * v22 * v27
          + v14 * v16 * sqrt(1 - v10 + v17) * v22 * v27 + v10 * v13 * v22 * v23 * v27 * v30
          + v11 * v13 * v22 * v23 * v27 * v30 + v11 * v12 * v20 * v21 * v26 * v31
          + v11 * v13 * v20 * v21 * v26 * v31 + v10 * v12 * v20 * v21 * v26 * v30 * v31
          + v10 * v13 * v20 * v21 * v26 * v30 * v31 + v10 * v12 * v22 * v23 * v27 * v30 * v31
          + v11 * v12 * v22 * v23 * v27 * v30 * v31 + v12 * v16 * v20 * v26 * v31 * v7
          + v13 * v16 * v20 * v26 * v31 * v7
          + sqrt(1 - v13 + v18) * v19 * v20 * v25 * v26 * v30 * v6 * v7
          + v18 * sqrt(1 - v12 + v19) * v20 * v25 * v26 * v30 * v6 * v7
          + v12 * v20 * v25 * v26 * v31 * v6 * v7 + v13 * v20 * v25 * v26 * v31 * v6 * v7
          + v12 * sqrt(1 - v11 + v16) * v17 * v22 * v27 * v9
          + v12 * v16 * sqrt(1 - v10 + v17) * v22 * v27 * v9 + v10 * v22 * v24 * v27 * v30 * v8 * v9
          + v11 * v22 * v24 * v27 * v30 * v8 * v9
          + sqrt(1 - v11 + v16) * v17 * v22 * v24 * v27 * v31 * v8 * v9
          + v16 * sqrt(1 - v10 + v17) * v22 * v24 * v27 * v31 * v8 * v9)
         / ((v10 + v11) * (v12 + v13) * sqrt(v28) * sqrt(v29));
}