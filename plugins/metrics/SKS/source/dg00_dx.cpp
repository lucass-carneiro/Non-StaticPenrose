#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::dllgSKS_00_dx(double t, double x, double y, double z) const noexcept
    -> double {

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
  const double v14 = dZ_dx(2, t, x, y, z);
  const double v15 = dZ_dx(1, t, x, y, z);
  const double v16 = dY_dx(2, t, x, y, z);
  const double v17 = dY_dx(1, t, x, y, z);
  const double v18 = dX_dx(2, t, x, y, z);
  const double v19 = dX_dx(1, t, x, y, z);
  const double v20 = l2_KS(a2, v2, v1, v0);
  const double v21 = l1_KS(a2, v2, v1, v0);
  const double v22 = l2_KS(a1, v5, v4, v3);
  const double v23 = l1_KS(a1, v5, v4, v3);
  const double v24 = H_KS(M2, a2, v2, v1, v0);
  const double v25 = H_KS(M1, a1, v5, v4, v3);
  const double v26 = pow(v20, 2);
  const double v27 = pow(v21, 2);
  const double v28 = pow(v22, 2);
  const double v29 = pow(v23, 2);
  const double v30 = dl2_KS_dZ(a2, v2, v1, v0);
  const double v31 = dl1_KS_dZ(a2, v2, v1, v0);
  const double v32 = dl2_KS_dZ(a1, v5, v4, v3);
  const double v33 = dl1_KS_dZ(a1, v5, v4, v3);
  const double v34 = dl2_KS_dY(a2, v2, v1, v0);
  const double v35 = dl1_KS_dY(a2, v2, v1, v0);
  const double v36 = dl2_KS_dY(a1, v5, v4, v3);
  const double v37 = dl1_KS_dY(a1, v5, v4, v3);
  const double v38 = dl2_KS_dX(a2, v2, v1, v0);
  const double v39 = dl1_KS_dX(a2, v2, v1, v0);
  const double v40 = dl2_KS_dX(a1, v5, v4, v3);
  const double v41 = dl1_KS_dX(a1, v5, v4, v3);
  const double v42 = dH_KS_dZ(M2, a2, v2, v1, v0);
  const double v43 = dH_KS_dZ(M1, a1, v5, v4, v3);
  const double v44 = dH_KS_dY(M2, a2, v2, v1, v0);
  const double v45 = dH_KS_dY(M1, a1, v5, v4, v3);
  const double v46 = dH_KS_dX(M2, a2, v2, v1, v0);
  const double v47 = dH_KS_dX(M1, a1, v5, v4, v3);

  return (2 * v10 * v14 * v20 * v24 * v30 - 2 * v10 * v11 * v14 * v20 * v24 * v30
          - 2 * v10 * v12 * v14 * v20 * v24 * v30 + 2 * v13 * v14 * v21 * v24 * v31
          - 2 * v11 * v13 * v14 * v21 * v24 * v31 - 2 * v12 * v13 * v14 * v21 * v24 * v31
          + 2 * v11 * v15 * v22 * v25 * v32 - 2 * v10 * v11 * v15 * v22 * v25 * v32
          - 2 * v11 * v13 * v15 * v22 * v25 * v32 + 2 * v12 * v15 * v23 * v25 * v33
          - 2 * v10 * v12 * v15 * v23 * v25 * v33 - 2 * v12 * v13 * v15 * v23 * v25 * v33
          + 2 * v10 * v16 * v20 * v24 * v34 - 2 * v10 * v11 * v16 * v20 * v24 * v34
          - 2 * v10 * v12 * v16 * v20 * v24 * v34 + 2 * v13 * v16 * v21 * v24 * v35
          - 2 * v11 * v13 * v16 * v21 * v24 * v35 - 2 * v12 * v13 * v16 * v21 * v24 * v35
          + 2 * v11 * v17 * v22 * v25 * v36 - 2 * v10 * v11 * v17 * v22 * v25 * v36
          - 2 * v11 * v13 * v17 * v22 * v25 * v36 + 2 * v12 * v17 * v23 * v25 * v37
          - 2 * v10 * v12 * v17 * v23 * v25 * v37 - 2 * v12 * v13 * v17 * v23 * v25 * v37
          + 2 * v10 * v18 * v20 * v24 * v38 - 2 * v10 * v11 * v18 * v20 * v24 * v38
          - 2 * v10 * v12 * v18 * v20 * v24 * v38 + 2 * v13 * v18 * v21 * v24 * v39
          - 2 * v11 * v13 * v18 * v21 * v24 * v39 - 2 * v12 * v13 * v18 * v21 * v24 * v39
          + 2 * v11 * v19 * v22 * v25 * v40 - 2 * v10 * v11 * v19 * v22 * v25 * v40
          - 2 * v11 * v13 * v19 * v22 * v25 * v40 + 2 * v12 * v19 * v23 * v25 * v41
          - 2 * v10 * v12 * v19 * v23 * v25 * v41 - 2 * v12 * v13 * v19 * v23 * v25 * v41
          + v14 * v42 - v11 * v14 * v42 - v12 * v14 * v42 + v10 * v14 * v26 * v42
          - v10 * v11 * v14 * v26 * v42 - v10 * v12 * v14 * v26 * v42 + v13 * v14 * v27 * v42
          - v11 * v13 * v14 * v27 * v42 - v12 * v13 * v14 * v27 * v42 + v15 * v43 - v10 * v15 * v43
          - v13 * v15 * v43 + v11 * v15 * v28 * v43 - v10 * v11 * v15 * v28 * v43
          - v11 * v13 * v15 * v28 * v43 + v12 * v15 * v29 * v43 - v10 * v12 * v15 * v29 * v43
          - v12 * v13 * v15 * v29 * v43 + v16 * v44 - v11 * v16 * v44 - v12 * v16 * v44
          + v10 * v16 * v26 * v44 - v10 * v11 * v16 * v26 * v44 - v10 * v12 * v16 * v26 * v44
          + v13 * v16 * v27 * v44 - v11 * v13 * v16 * v27 * v44 - v12 * v13 * v16 * v27 * v44
          + v17 * v45 - v10 * v17 * v45 - v13 * v17 * v45 + v11 * v17 * v28 * v45
          - v10 * v11 * v17 * v28 * v45 - v11 * v13 * v17 * v28 * v45 + v12 * v17 * v29 * v45
          - v10 * v12 * v17 * v29 * v45 - v12 * v13 * v17 * v29 * v45 + v18 * v46 - v11 * v18 * v46
          - v12 * v18 * v46 + v10 * v18 * v26 * v46 - v10 * v11 * v18 * v26 * v46
          - v10 * v12 * v18 * v26 * v46 + v13 * v18 * v27 * v46 - v11 * v13 * v18 * v27 * v46
          - v12 * v13 * v18 * v27 * v46 + v19 * v47 - v10 * v19 * v47 - v13 * v19 * v47
          + v11 * v19 * v28 * v47 - v10 * v11 * v19 * v28 * v47 - v11 * v13 * v19 * v28 * v47
          + v12 * v19 * v29 * v47 - v10 * v12 * v19 * v29 * v47 - v12 * v13 * v19 * v29 * v47
          - 2 * v14 * v24 * v30 * v6 + 2 * v11 * v14 * v24 * v30 * v6
          + 2 * v12 * v14 * v24 * v30 * v6 - 2 * v16 * v24 * v34 * v6
          + 2 * v11 * v16 * v24 * v34 * v6 + 2 * v12 * v16 * v24 * v34 * v6
          - 2 * v18 * v24 * v38 * v6 + 2 * v11 * v18 * v24 * v38 * v6
          + 2 * v12 * v18 * v24 * v38 * v6 - 2 * v14 * v20 * v42 * v6
          + 2 * v11 * v14 * v20 * v42 * v6 + 2 * v12 * v14 * v20 * v42 * v6
          - 2 * v16 * v20 * v44 * v6 + 2 * v11 * v16 * v20 * v44 * v6
          + 2 * v12 * v16 * v20 * v44 * v6 - 2 * v18 * v20 * v46 * v6
          + 2 * v11 * v18 * v20 * v46 * v6 + 2 * v12 * v18 * v20 * v46 * v6
          - 2 * v15 * v25 * v32 * v7 + 2 * v10 * v15 * v25 * v32 * v7
          + 2 * v13 * v15 * v25 * v32 * v7 - 2 * v17 * v25 * v36 * v7
          + 2 * v10 * v17 * v25 * v36 * v7 + 2 * v13 * v17 * v25 * v36 * v7
          - 2 * v19 * v25 * v40 * v7 + 2 * v10 * v19 * v25 * v40 * v7
          + 2 * v13 * v19 * v25 * v40 * v7 - 2 * v15 * v22 * v43 * v7
          + 2 * v10 * v15 * v22 * v43 * v7 + 2 * v13 * v15 * v22 * v43 * v7
          - 2 * v17 * v22 * v45 * v7 + 2 * v10 * v17 * v22 * v45 * v7
          + 2 * v13 * v17 * v22 * v45 * v7 - 2 * v19 * v22 * v47 * v7
          + 2 * v10 * v19 * v22 * v47 * v7 + 2 * v13 * v19 * v22 * v47 * v7
          - 2 * v15 * v25 * v33 * v8 + 2 * v10 * v15 * v25 * v33 * v8
          + 2 * v13 * v15 * v25 * v33 * v8 - 2 * v17 * v25 * v37 * v8
          + 2 * v10 * v17 * v25 * v37 * v8 + 2 * v13 * v17 * v25 * v37 * v8
          - 2 * v19 * v25 * v41 * v8 + 2 * v10 * v19 * v25 * v41 * v8
          + 2 * v13 * v19 * v25 * v41 * v8 - 2 * v15 * v23 * v43 * v8
          + 2 * v10 * v15 * v23 * v43 * v8 + 2 * v13 * v15 * v23 * v43 * v8
          - 2 * v17 * v23 * v45 * v8 + 2 * v10 * v17 * v23 * v45 * v8
          + 2 * v13 * v17 * v23 * v45 * v8 - 2 * v19 * v23 * v47 * v8
          + 2 * v10 * v19 * v23 * v47 * v8 + 2 * v13 * v19 * v23 * v47 * v8
          + 2 * v15 * v23 * v25 * v32 * v7 * v8 - 2 * v10 * v15 * v23 * v25 * v32 * v7 * v8
          - 2 * v13 * v15 * v23 * v25 * v32 * v7 * v8 + 2 * v15 * v22 * v25 * v33 * v7 * v8
          - 2 * v10 * v15 * v22 * v25 * v33 * v7 * v8 - 2 * v13 * v15 * v22 * v25 * v33 * v7 * v8
          + 2 * v17 * v23 * v25 * v36 * v7 * v8 - 2 * v10 * v17 * v23 * v25 * v36 * v7 * v8
          - 2 * v13 * v17 * v23 * v25 * v36 * v7 * v8 + 2 * v17 * v22 * v25 * v37 * v7 * v8
          - 2 * v10 * v17 * v22 * v25 * v37 * v7 * v8 - 2 * v13 * v17 * v22 * v25 * v37 * v7 * v8
          + 2 * v19 * v23 * v25 * v40 * v7 * v8 - 2 * v10 * v19 * v23 * v25 * v40 * v7 * v8
          - 2 * v13 * v19 * v23 * v25 * v40 * v7 * v8 + 2 * v19 * v22 * v25 * v41 * v7 * v8
          - 2 * v10 * v19 * v22 * v25 * v41 * v7 * v8 - 2 * v13 * v19 * v22 * v25 * v41 * v7 * v8
          + 2 * v15 * v22 * v23 * v43 * v7 * v8 - 2 * v10 * v15 * v22 * v23 * v43 * v7 * v8
          - 2 * v13 * v15 * v22 * v23 * v43 * v7 * v8 + 2 * v17 * v22 * v23 * v45 * v7 * v8
          - 2 * v10 * v17 * v22 * v23 * v45 * v7 * v8 - 2 * v13 * v17 * v22 * v23 * v45 * v7 * v8
          + 2 * v19 * v22 * v23 * v47 * v7 * v8 - 2 * v10 * v19 * v22 * v23 * v47 * v7 * v8
          - 2 * v13 * v19 * v22 * v23 * v47 * v7 * v8 - 2 * v14 * v24 * v31 * v9
          + 2 * v11 * v14 * v24 * v31 * v9 + 2 * v12 * v14 * v24 * v31 * v9
          - 2 * v16 * v24 * v35 * v9 + 2 * v11 * v16 * v24 * v35 * v9
          + 2 * v12 * v16 * v24 * v35 * v9 - 2 * v18 * v24 * v39 * v9
          + 2 * v11 * v18 * v24 * v39 * v9 + 2 * v12 * v18 * v24 * v39 * v9
          - 2 * v14 * v21 * v42 * v9 + 2 * v11 * v14 * v21 * v42 * v9
          + 2 * v12 * v14 * v21 * v42 * v9 - 2 * v16 * v21 * v44 * v9
          + 2 * v11 * v16 * v21 * v44 * v9 + 2 * v12 * v16 * v21 * v44 * v9
          - 2 * v18 * v21 * v46 * v9 + 2 * v11 * v18 * v21 * v46 * v9
          + 2 * v12 * v18 * v21 * v46 * v9 + 2 * v14 * v21 * v24 * v30 * v6 * v9
          - 2 * v11 * v14 * v21 * v24 * v30 * v6 * v9 - 2 * v12 * v14 * v21 * v24 * v30 * v6 * v9
          + 2 * v14 * v20 * v24 * v31 * v6 * v9 - 2 * v11 * v14 * v20 * v24 * v31 * v6 * v9
          - 2 * v12 * v14 * v20 * v24 * v31 * v6 * v9 + 2 * v16 * v21 * v24 * v34 * v6 * v9
          - 2 * v11 * v16 * v21 * v24 * v34 * v6 * v9 - 2 * v12 * v16 * v21 * v24 * v34 * v6 * v9
          + 2 * v16 * v20 * v24 * v35 * v6 * v9 - 2 * v11 * v16 * v20 * v24 * v35 * v6 * v9
          - 2 * v12 * v16 * v20 * v24 * v35 * v6 * v9 + 2 * v18 * v21 * v24 * v38 * v6 * v9
          - 2 * v11 * v18 * v21 * v24 * v38 * v6 * v9 - 2 * v12 * v18 * v21 * v24 * v38 * v6 * v9
          + 2 * v18 * v20 * v24 * v39 * v6 * v9 - 2 * v11 * v18 * v20 * v24 * v39 * v6 * v9
          - 2 * v12 * v18 * v20 * v24 * v39 * v6 * v9 + 2 * v14 * v20 * v21 * v42 * v6 * v9
          - 2 * v11 * v14 * v20 * v21 * v42 * v6 * v9 - 2 * v12 * v14 * v20 * v21 * v42 * v6 * v9
          + 2 * v16 * v20 * v21 * v44 * v6 * v9 - 2 * v11 * v16 * v20 * v21 * v44 * v6 * v9
          - 2 * v12 * v16 * v20 * v21 * v44 * v6 * v9 + 2 * v18 * v20 * v21 * v46 * v6 * v9
          - 2 * v11 * v18 * v20 * v21 * v46 * v6 * v9 - 2 * v12 * v18 * v20 * v21 * v46 * v6 * v9)
         / ((-1 + v11 + v12) * (-1 + v10 + v13));
}