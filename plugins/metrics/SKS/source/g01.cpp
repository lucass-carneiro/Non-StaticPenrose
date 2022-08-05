#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_01(double t, double x, double y, double z) const noexcept -> double {

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
  const double v12 = pow(v6, 3);
  const double v13 = pow(v8, 4);
  const double v14 = pow(v8, 2);
  const double v15 = pow(v9, 2);
  const double v16 = pow(v9, 4);
  const double v17 = pow(v6, 4);
  const double v18 = pow(v8, 3);
  const double v19 = pow(v9, 3);
  const double v20 = pow(v7, 3);
  const double v21 = pow(v7, 4);
  const double v22 = -v10;
  const double v23 = -v11;
  const double v24 = -v14;
  const double v25 = -v15;
  const double v26 = l2_KS(a2, v2, v1, v0);
  const double v27 = l1_KS(a2, v2, v1, v0);
  const double v28 = l2_KS(a1, v5, v4, v3);
  const double v29 = l1_KS(a1, v5, v4, v3);
  const double v30 = H_KS(M2, a2, v2, v1, v0);
  const double v31 = H_KS(M1, a1, v5, v4, v3);
  const double v32 = pow(v26, 2);
  const double v33 = pow(v27, 2);
  const double v34 = pow(v28, 2);
  const double v35 = pow(v29, 2);
  const double v36 = sqrt(1 + v22 + v23);
  const double v37 = sqrt(1 + v24 + v25);

  return (v13 * v20 * v30 + 2 * v14 * v15 * v20 * v30 + v16 * v20 * v30 + v20 * v24 * v30
          + v20 * v25 * v30 + v11 * v14 * v27 * v30 + v10 * v11 * v14 * v27 * v30
          + v11 * v15 * v27 * v30 + v10 * v11 * v15 * v27 * v30 - 2 * v11 * v14 * v15 * v27 * v30
          - 2 * v10 * v11 * v14 * v15 * v27 * v30 - v13 * v21 * v27 * v30 + v14 * v21 * v27 * v30
          + v15 * v21 * v27 * v30 - 2 * v14 * v15 * v21 * v27 * v30 - v16 * v21 * v27 * v30
          + v11 * v13 * v22 * v27 * v30 + v11 * v16 * v22 * v27 * v30 + v13 * v23 * v27 * v30
          + v16 * v23 * v27 * v30 + v13 * v22 * sqrt(1 - v10 + v23) * v27 * v30
          + v16 * v22 * sqrt(1 - v10 + v23) * v27 * v30 + 2 * v10 * v11 * v19 * v31
          + v17 * v19 * v31 + v19 * v21 * v31 + v19 * v22 * v31 + v19 * v23 * v31
          + v10 * v15 * v29 * v31 + v11 * v15 * v29 * v31 - 2 * v10 * v11 * v15 * v29 * v31
          + v10 * v14 * v15 * v29 * v31 + v11 * v14 * v15 * v29 * v31
          - 2 * v10 * v11 * v14 * v15 * v29 * v31 + v10 * v16 * v29 * v31 + v11 * v16 * v29 * v31
          - 2 * v10 * v11 * v16 * v29 * v31 - v16 * v17 * v29 * v31 - v16 * v21 * v29 * v31
          + v15 * v17 * v24 * v29 * v31 + v15 * v21 * v24 * v29 * v31 + v17 * v25 * v29 * v31
          + v21 * v25 * v29 * v31 + v17 * v24 * sqrt(1 - v14 + v25) * v29 * v31
          + v21 * v24 * sqrt(1 - v14 + v25) * v29 * v31 + v13 * v20 * v30 * v33
          + 2 * v14 * v15 * v20 * v30 * v33 + v16 * v20 * v30 * v33 + v20 * v24 * v30 * v33
          + v20 * v25 * v30 * v33 + 2 * v10 * v11 * v19 * v31 * v35 + v17 * v19 * v31 * v35
          + v19 * v21 * v31 * v35 + v19 * v22 * v31 * v35 + v19 * v23 * v31 * v35
          + v10 * v14 * v27 * v30 * v36 + v10 * v15 * v27 * v30 * v36
          - 2 * v10 * v14 * v15 * v27 * v30 * v36 + v12 * v13 * v26 * v27 * v30 * v36
          + 2 * v12 * v14 * v15 * v26 * v27 * v30 * v36 + v12 * v16 * v26 * v27 * v30 * v36
          + v12 * v24 * v26 * v27 * v30 * v36 + v12 * v25 * v26 * v27 * v30 * v36
          + v10 * v14 * v29 * v31 * v37 + v11 * v14 * v29 * v31 * v37
          - 2 * v10 * v11 * v14 * v29 * v31 * v37 + 2 * v10 * v11 * v18 * v28 * v29 * v31 * v37
          + v17 * v18 * v28 * v29 * v31 * v37 + v18 * v21 * v28 * v29 * v31 * v37
          + v18 * v22 * v28 * v29 * v31 * v37 + v18 * v23 * v28 * v29 * v31 * v37
          - v13 * v20 * v26 * v30 * v6 + v14 * v20 * v26 * v30 * v6 + v15 * v20 * v26 * v30 * v6
          - 2 * v14 * v15 * v20 * v26 * v30 * v6 - v16 * v20 * v26 * v30 * v6
          + 2 * v11 * v13 * v26 * v27 * v30 * v6 - 2 * v11 * v14 * v26 * v27 * v30 * v6
          - 2 * v11 * v15 * v26 * v27 * v30 * v6 + 4 * v11 * v14 * v15 * v26 * v27 * v30 * v6
          + 2 * v11 * v16 * v26 * v27 * v30 * v6
          + v13 * sqrt(1 - v11 + v22) * v23 * v26 * v27 * v30 * v6
          + v16 * sqrt(1 - v11 + v22) * v23 * v26 * v27 * v30 * v6
          + v11 * v14 * v26 * v27 * v30 * v36 * v6 + v11 * v15 * v26 * v27 * v30 * v36 * v6
          - 2 * v11 * v14 * v15 * v26 * v27 * v30 * v36 * v6 + v10 * v13 * v30 * v7
          + 2 * v10 * v14 * v15 * v30 * v7 + v10 * v16 * v30 * v7 + v14 * v22 * v30 * v7
          + v15 * v22 * v30 * v7 - v12 * v13 * v26 * v30 * v7 + v12 * v14 * v26 * v30 * v7
          + v12 * v15 * v26 * v30 * v7 - 2 * v12 * v14 * v15 * v26 * v30 * v7
          - v12 * v16 * v26 * v30 * v7 + v10 * v13 * v30 * v32 * v7
          + 2 * v10 * v14 * v15 * v30 * v32 * v7 + v10 * v16 * v30 * v32 * v7
          + v14 * v22 * v30 * v32 * v7 + v15 * v22 * v30 * v32 * v7
          + v13 * v22 * sqrt(1 - v10 + v23) * v30 * v32 * v7
          + v16 * v22 * sqrt(1 - v10 + v23) * v30 * v32 * v7
          + v14 * v22 * sqrt(1 - v10 + v23) * v30 * v33 * v7
          + v15 * v22 * sqrt(1 - v10 + v23) * v30 * v33 * v7 + v10 * v14 * v30 * v32 * v36 * v7
          + v10 * v15 * v30 * v32 * v36 * v7 - 2 * v10 * v14 * v15 * v30 * v32 * v36 * v7
          + v10 * v13 * v30 * v33 * v36 * v7 + 2 * v10 * v14 * v15 * v30 * v33 * v36 * v7
          + v10 * v16 * v30 * v33 * v36 * v7 - v13 * v26 * v30 * v6 * v7 + v14 * v26 * v30 * v6 * v7
          + v15 * v26 * v30 * v6 * v7 - 2 * v14 * v15 * v26 * v30 * v6 * v7
          - v16 * v26 * v30 * v6 * v7 + v13 * v26 * v30 * v36 * v6 * v7
          + 2 * v14 * v15 * v26 * v30 * v36 * v6 * v7 + v16 * v26 * v30 * v36 * v6 * v7
          + v24 * v26 * v30 * v36 * v6 * v7 + v25 * v26 * v30 * v36 * v6 * v7
          + v10 * v19 * v28 * v31 * v8 + v11 * v19 * v28 * v31 * v8
          - 2 * v10 * v11 * v19 * v28 * v31 * v8 - v17 * v19 * v28 * v31 * v8
          - v19 * v21 * v28 * v31 * v8 - 2 * v10 * v15 * v28 * v29 * v31 * v8
          - 2 * v11 * v15 * v28 * v29 * v31 * v8 + 4 * v10 * v11 * v15 * v28 * v29 * v31 * v8
          + 2 * v15 * v17 * v28 * v29 * v31 * v8 + 2 * v15 * v21 * v28 * v29 * v31 * v8
          + v17 * sqrt(1 - v15 + v24) * v25 * v28 * v29 * v31 * v8
          + v21 * sqrt(1 - v15 + v24) * v25 * v28 * v29 * v31 * v8
          + v10 * v15 * v28 * v29 * v31 * v37 * v8 + v11 * v15 * v28 * v29 * v31 * v37 * v8
          - 2 * v10 * v11 * v15 * v28 * v29 * v31 * v37 * v8 + 2 * v10 * v11 * v14 * v31 * v9
          + v14 * v17 * v31 * v9 + v14 * v21 * v31 * v9 + v14 * v22 * v31 * v9
          + v14 * v23 * v31 * v9 + v10 * v18 * v28 * v31 * v9 + v11 * v18 * v28 * v31 * v9
          - 2 * v10 * v11 * v18 * v28 * v31 * v9 - v17 * v18 * v28 * v31 * v9
          - v18 * v21 * v28 * v31 * v9 + 2 * v10 * v11 * v14 * v31 * v34 * v9
          + v14 * v17 * v31 * v34 * v9 + v14 * v21 * v31 * v34 * v9 + v14 * v22 * v31 * v34 * v9
          + v14 * v23 * v31 * v34 * v9 + v17 * v24 * sqrt(1 - v14 + v25) * v31 * v34 * v9
          + v21 * v24 * sqrt(1 - v14 + v25) * v31 * v34 * v9 + v10 * v14 * v31 * v34 * v37 * v9
          + v11 * v14 * v31 * v34 * v37 * v9 - 2 * v10 * v11 * v14 * v31 * v34 * v37 * v9
          + 2 * v10 * v11 * v14 * v31 * v35 * v37 * v9 + v14 * v17 * v31 * v35 * v37 * v9
          + v14 * v21 * v31 * v35 * v37 * v9 + v14 * v22 * v31 * v35 * v37 * v9
          + v14 * v23 * v31 * v35 * v37 * v9 + v10 * v28 * v31 * v8 * v9 + v11 * v28 * v31 * v8 * v9
          - 2 * v10 * v11 * v28 * v31 * v8 * v9 - v17 * v28 * v31 * v8 * v9
          - v21 * v28 * v31 * v8 * v9 + 2 * v10 * v11 * v28 * v31 * v37 * v8 * v9
          + v17 * v28 * v31 * v37 * v8 * v9 + v21 * v28 * v31 * v37 * v8 * v9
          + v22 * v28 * v31 * v37 * v8 * v9 + v23 * v28 * v31 * v37 * v8 * v9)
         / ((-1 + v10 + v11) * (v10 + v11) * (-1 + v14 + v15) * (v14 + v15));
}