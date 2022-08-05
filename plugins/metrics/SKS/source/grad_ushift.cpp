#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::grad_ushift(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {

  using std::pow;

  metric_server::spatial_matrix gradushift{};

  const double v0 = llgSKS_23(t, x, y, z);
  const double v1 = llgSKS_22(t, x, y, z);
  const double v2 = llgSKS_12(t, x, y, z);
  const double v3 = llgSKS_11(t, x, y, z);
  const double v4 = llgSKS_01(t, x, y, z);
  const double v5 = llgSKS_02(t, x, y, z);
  const double v6 = llgSKS_13(t, x, y, z);
  const double v7 = llgSKS_03(t, x, y, z);
  const double v8 = llgSKS_33(t, x, y, z);
  const double v9 = pow(v3, 2);
  const double v10 = pow(v2, 3);
  const double v11 = pow(v2, 2);
  const double v12 = pow(v1, 2);
  const double v13 = pow(v2, 4);
  const double v14 = pow(v0, 2);
  const double v15 = pow(v6, 2);
  const double v16 = pow(v6, 3);
  const double v17 = pow(v0, 3);
  const double v18 = pow(v8, 2);
  const double v19 = pow(v6, 4);
  const double v20 = pow(v0, 4);
  const double v21 = dllgSKS_33_dz(t, x, y, z);
  const double v22 = dllgSKS_23_dz(t, x, y, z);
  const double v23 = dllgSKS_22_dz(t, x, y, z);
  const double v24 = dllgSKS_13_dz(t, x, y, z);
  const double v25 = dllgSKS_12_dz(t, x, y, z);
  const double v26 = dllgSKS_11_dz(t, x, y, z);
  const double v27 = dllgSKS_03_dz(t, x, y, z);
  const double v28 = dllgSKS_02_dz(t, x, y, z);
  const double v29 = dllgSKS_01_dz(t, x, y, z);
  const double v30 = dllgSKS_33_dy(t, x, y, z);
  const double v31 = dllgSKS_23_dy(t, x, y, z);
  const double v32 = dllgSKS_22_dy(t, x, y, z);
  const double v33 = dllgSKS_13_dy(t, x, y, z);
  const double v34 = dllgSKS_12_dy(t, x, y, z);
  const double v35 = dllgSKS_11_dy(t, x, y, z);
  const double v36 = dllgSKS_03_dy(t, x, y, z);
  const double v37 = dllgSKS_02_dy(t, x, y, z);
  const double v38 = dllgSKS_01_dy(t, x, y, z);
  const double v39 = dllgSKS_33_dx(t, x, y, z);
  const double v40 = dllgSKS_23_dx(t, x, y, z);
  const double v41 = dllgSKS_22_dx(t, x, y, z);
  const double v42 = dllgSKS_13_dx(t, x, y, z);
  const double v43 = dllgSKS_12_dx(t, x, y, z);
  const double v44 = dllgSKS_11_dx(t, x, y, z);
  const double v45 = dllgSKS_03_dx(t, x, y, z);
  const double v46 = dllgSKS_02_dx(t, x, y, z);
  const double v47 = dllgSKS_01_dx(t, x, y, z);
  const double v48 = v11 * v8;
  const double v49 = v14 * v3;
  const double v50 = v1 * v15;
  const double v51 = pow(v48 + v49 + v50 - 2 * v0 * v2 * v6 - v1 * v3 * v8, -2);

  gradushift[0][0]
      = v51
        * (-(v11 * v14 * v39 * v4) - v12 * v15 * v39 * v4 - v14 * v15 * v4 * v41
           - v11 * v18 * v4 * v41 + 2 * v17 * v2 * v4 * v42 + 2 * v1 * v18 * v2 * v4 * v43
           - v12 * v18 * v4 * v44 - v20 * v4 * v44 + v12 * v16 * v45 - v17 * v2 * v3 * v45
           - v0 * v1 * v16 * v46 + v10 * v18 * v46 + 2 * v14 * v15 * v2 * v46
           - v1 * v18 * v2 * v3 * v46 - v1 * v11 * v18 * v47 + v12 * v18 * v3 * v47 + v20 * v3 * v47
           + 2 * v0 * v4 * v40 * v48 + v14 * v47 * v48 - v1 * v16 * v40 * v5 + v0 * v16 * v41 * v5
           + v18 * v2 * v3 * v41 * v5 - v17 * v3 * v42 * v5 - 2 * v14 * v15 * v43 * v5
           - v11 * v18 * v43 * v5 - v1 * v18 * v3 * v43 * v5 + v1 * v18 * v2 * v44 * v5
           + v0 * v42 * v48 * v5 + v2 * v39 * v49 * v5 + 2 * v0 * v4 * v40 * v50
           - 3 * v0 * v2 * v45 * v50 + v14 * v47 * v50 + v2 * v39 * v5 * v50 + v0 * v42 * v5 * v50
           + 2 * v0 * v1 * v2 * v39 * v4 * v6 - 2 * v14 * v2 * v4 * v40 * v6
           - 2 * v1 * v14 * v4 * v42 * v6 + 2 * v17 * v4 * v43 * v6 + 2 * v11 * v14 * v45 * v6
           - v17 * v3 * v46 * v6 - 2 * v17 * v2 * v47 * v6 + v1 * v45 * v48 * v6
           - 3 * v0 * v46 * v48 * v6 + v1 * v45 * v49 * v6 - v0 * v11 * v39 * v5 * v6
           - v0 * v1 * v3 * v39 * v5 * v6 + v17 * v44 * v5 * v6 + v40 * v48 * v5 * v6
           + v40 * v49 * v5 * v6 + v0 * v10 * v39 * v7 - v0 * v1 * v2 * v3 * v39 * v7
           - v0 * v15 * v2 * v41 * v7 - 2 * v11 * v14 * v42 * v7 - v12 * v15 * v42 * v7
           - v17 * v3 * v43 * v7 + v17 * v2 * v44 * v7 + v1 * v42 * v48 * v7 + v0 * v43 * v48 * v7
           + v2 * v40 * v49 * v7 + v1 * v42 * v49 * v7 + v2 * v40 * v50 * v7 + v0 * v43 * v50 * v7
           - v1 * v11 * v39 * v6 * v7 + v12 * v3 * v39 * v6 * v7 - 2 * v0 * v1 * v3 * v40 * v6 * v7
           + 2 * v0 * v1 * v2 * v42 * v6 * v7 - v1 * v14 * v44 * v6 * v7 + v41 * v48 * v6 * v7
           + v41 * v49 * v6 * v7 - 2 * v0 * v1 * v2 * v4 * v42 * v8 - 2 * v14 * v2 * v4 * v43 * v8
           + 2 * v1 * v14 * v4 * v44 * v8 - v0 * v10 * v45 * v8 + v0 * v1 * v2 * v3 * v45 * v8
           - v12 * v15 * v47 * v8 + v2 * v46 * v49 * v8 - 2 * v1 * v47 * v49 * v8
           - 2 * v0 * v2 * v3 * v40 * v5 * v8 - v15 * v2 * v41 * v5 * v8
           + v0 * v1 * v3 * v42 * v5 * v8 - v14 * v2 * v44 * v5 * v8 + v43 * v49 * v5 * v8
           + v2 * v46 * v50 * v8 + v43 * v5 * v50 * v8 - 2 * v1 * v2 * v4 * v40 * v6 * v8
           + 2 * v0 * v2 * v4 * v41 * v6 * v8 + 2 * v12 * v4 * v42 * v6 * v8
           - 2 * v0 * v1 * v4 * v43 * v6 * v8 - v12 * v3 * v45 * v6 * v8
           + v0 * v1 * v3 * v46 * v6 * v8 + 2 * v0 * v1 * v2 * v47 * v6 * v8
           + v1 * v3 * v40 * v5 * v6 * v8 - v0 * v3 * v41 * v5 * v6 * v8
           - 2 * v1 * v2 * v42 * v5 * v6 * v8 + 2 * v0 * v2 * v43 * v5 * v6 * v8
           - v0 * v1 * v44 * v5 * v6 * v8 - v10 * v40 * v7 * v8 + v1 * v2 * v3 * v40 * v7 * v8
           - v0 * v2 * v3 * v41 * v7 * v8 - v12 * v3 * v42 * v7 * v8 + v0 * v1 * v3 * v43 * v7 * v8
           - v0 * v1 * v2 * v44 * v7 * v8 - 2 * v1 * v2 * v43 * v6 * v7 * v8
           + v12 * v44 * v6 * v7 * v8);

  gradushift[0][1]
      = v51
        * (-(v1 * v16 * v4 * v40) + v0 * v16 * v4 * v41 + v18 * v2 * v3 * v4 * v41
           - v17 * v3 * v4 * v42 - 2 * v14 * v15 * v4 * v43 - v11 * v18 * v4 * v43
           - v1 * v18 * v3 * v4 * v43 + v1 * v18 * v2 * v4 * v44 + 2 * v0 * v11 * v15 * v45
           - v1 * v16 * v2 * v45 + v1 * v19 * v46 - 2 * v0 * v16 * v2 * v46 - v11 * v18 * v3 * v46
           - v0 * v1 * v16 * v47 + v10 * v18 * v47 + 2 * v14 * v15 * v2 * v47
           - v1 * v18 * v2 * v3 * v47 + v0 * v4 * v42 * v48 + v0 * v3 * v45 * v48 + v15 * v46 * v48
           + v2 * v39 * v4 * v49 + v15 * v46 * v49 - v11 * v15 * v39 * v5 + 2 * v16 * v2 * v40 * v5
           - 2 * v0 * v15 * v3 * v40 * v5 - v19 * v41 * v5 - 2 * v0 * v15 * v2 * v42 * v5
           + 2 * v0 * v16 * v43 * v5 + 2 * v18 * v2 * v3 * v43 * v5 - v14 * v15 * v44 * v5
           - v11 * v18 * v44 * v5 + v2 * v39 * v4 * v50 + v0 * v4 * v42 * v50 + v0 * v3 * v45 * v50
           - v0 * v11 * v39 * v4 * v6 - v0 * v1 * v3 * v39 * v4 * v6 + v17 * v4 * v44 * v6
           - v17 * v3 * v47 * v6 + v4 * v40 * v48 * v6 - 3 * v0 * v47 * v48 * v6
           + v4 * v40 * v49 * v6 - 3 * v2 * v45 * v49 * v6 + 2 * v0 * v2 * v3 * v39 * v5 * v6
           + 2 * v42 * v48 * v5 * v6 + 2 * v42 * v49 * v5 * v6 - v0 * v11 * v3 * v39 * v7
           - 2 * v11 * v15 * v40 * v7 + v16 * v2 * v41 * v7 - v0 * v15 * v3 * v41 * v7
           - v1 * v16 * v43 * v7 + v3 * v40 * v48 * v7 + v0 * v44 * v48 * v7 + v2 * v42 * v49 * v7
           + v3 * v40 * v50 * v7 + v2 * v42 * v50 * v7 + v0 * v44 * v50 * v7 + v10 * v39 * v6 * v7
           - v1 * v2 * v3 * v39 * v6 * v7 + 2 * v0 * v2 * v3 * v40 * v6 * v7
           - 2 * v0 * v1 * v3 * v42 * v6 * v7 - v14 * v2 * v44 * v6 * v7 + v43 * v48 * v6 * v7
           + v43 * v49 * v6 * v7 - 2 * v0 * v2 * v3 * v4 * v40 * v8 - v15 * v2 * v4 * v41 * v8
           + v0 * v1 * v3 * v4 * v42 * v8 - v14 * v2 * v4 * v44 * v8 + v4 * v43 * v49 * v8
           + v2 * v47 * v49 * v8 + 2 * v15 * v3 * v41 * v5 * v8 - 2 * v0 * v2 * v3 * v42 * v5 * v8
           - 2 * v15 * v2 * v43 * v5 * v8 + v4 * v43 * v50 * v8 - 2 * v3 * v46 * v50 * v8
           + v2 * v47 * v50 * v8 + v1 * v3 * v4 * v40 * v6 * v8 - v0 * v3 * v4 * v41 * v6 * v8
           - 2 * v1 * v2 * v4 * v42 * v6 * v8 + 2 * v0 * v2 * v4 * v43 * v6 * v8
           - v0 * v1 * v4 * v44 * v6 * v8 - v10 * v45 * v6 * v8 + v1 * v2 * v3 * v45 * v6 * v8
           + 2 * v0 * v2 * v3 * v46 * v6 * v8 + v0 * v1 * v3 * v47 * v6 * v8
           - 2 * v2 * v3 * v40 * v5 * v6 * v8 - 2 * v0 * v3 * v43 * v5 * v6 * v8
           + 2 * v0 * v2 * v44 * v5 * v6 * v8 - v10 * v42 * v7 * v8 + v1 * v2 * v3 * v42 * v7 * v8
           - 2 * v0 * v2 * v3 * v43 * v7 * v8 - v2 * v3 * v41 * v6 * v7 * v8
           + v1 * v3 * v43 * v6 * v7 * v8 - v1 * v2 * v44 * v6 * v7 * v8 + v17 * v45 * v9
           + v1 * v18 * v46 * v9 - v14 * v39 * v5 * v9 - v18 * v41 * v5 * v9
           + v0 * v1 * v39 * v7 * v9 - v14 * v40 * v7 * v9 - v0 * v1 * v45 * v8 * v9
           - v14 * v46 * v8 * v9 + 2 * v0 * v40 * v5 * v8 * v9 - v1 * v40 * v7 * v8 * v9
           + v0 * v41 * v7 * v8 * v9);

  gradushift[0][2]
      = v51
        * (v0 * v10 * v39 * v4 - v0 * v1 * v2 * v3 * v39 * v4 - v0 * v15 * v2 * v4 * v41
           - 2 * v11 * v14 * v4 * v42 - v12 * v15 * v4 * v42 - v17 * v3 * v4 * v43
           + v17 * v2 * v4 * v44 - v12 * v15 * v3 * v45 + 2 * v0 * v11 * v15 * v46
           - v1 * v16 * v2 * v46 + v12 * v16 * v47 - v17 * v2 * v3 * v47 + v1 * v4 * v42 * v48
           + v0 * v4 * v43 * v48 - 2 * v1 * v3 * v45 * v48 + v0 * v3 * v46 * v48
           + v2 * v4 * v40 * v49 + v1 * v4 * v42 * v49 + v11 * v45 * v49 - v0 * v11 * v3 * v39 * v5
           - 2 * v11 * v15 * v40 * v5 + v16 * v2 * v41 * v5 - v0 * v15 * v3 * v41 * v5
           - v1 * v16 * v43 * v5 + v3 * v40 * v48 * v5 + v0 * v44 * v48 * v5 + v2 * v42 * v49 * v5
           + v2 * v4 * v40 * v50 + v0 * v4 * v43 * v50 + v11 * v45 * v50 + v0 * v3 * v46 * v50
           - 3 * v0 * v2 * v47 * v50 + v3 * v40 * v5 * v50 + v2 * v42 * v5 * v50
           + v0 * v44 * v5 * v50 - v1 * v11 * v39 * v4 * v6 + v12 * v3 * v39 * v4 * v6
           - 2 * v0 * v1 * v3 * v4 * v40 * v6 + 2 * v0 * v1 * v2 * v4 * v42 * v6
           - v1 * v14 * v4 * v44 * v6 - 2 * v0 * v10 * v45 * v6 + 2 * v0 * v1 * v2 * v3 * v45 * v6
           + 2 * v11 * v14 * v47 * v6 + v4 * v41 * v48 * v6 + v1 * v47 * v48 * v6
           + v4 * v41 * v49 * v6 - 3 * v2 * v46 * v49 * v6 + v1 * v47 * v49 * v6
           + v10 * v39 * v5 * v6 - v1 * v2 * v3 * v39 * v5 * v6 + 2 * v0 * v2 * v3 * v40 * v5 * v6
           - 2 * v0 * v1 * v3 * v42 * v5 * v6 - v14 * v2 * v44 * v5 * v6 + v43 * v48 * v5 * v6
           + v43 * v49 * v5 * v6 - v13 * v39 * v7 + 2 * v1 * v11 * v3 * v39 * v7
           - 2 * v0 * v11 * v3 * v40 * v7 - v11 * v15 * v41 * v7 + 2 * v0 * v10 * v42 * v7
           - 2 * v0 * v1 * v2 * v3 * v42 * v7 - v11 * v14 * v44 * v7 - v12 * v15 * v44 * v7
           + 2 * v2 * v43 * v49 * v7 + 2 * v2 * v43 * v50 * v7 + 2 * v10 * v40 * v6 * v7
           - 2 * v1 * v2 * v3 * v40 * v6 * v7 + 2 * v0 * v2 * v3 * v41 * v6 * v7
           - 2 * v1 * v11 * v42 * v6 * v7 + 2 * v12 * v3 * v42 * v6 * v7
           - 2 * v0 * v11 * v43 * v6 * v7 - 2 * v0 * v1 * v3 * v43 * v6 * v7
           + 2 * v0 * v1 * v2 * v44 * v6 * v7 - v10 * v4 * v40 * v8 + v1 * v2 * v3 * v4 * v40 * v8
           - v0 * v2 * v3 * v4 * v41 * v8 - v12 * v3 * v4 * v42 * v8 + v0 * v1 * v3 * v4 * v43 * v8
           - v0 * v1 * v2 * v4 * v44 * v8 + v13 * v45 * v8 - v0 * v10 * v47 * v8
           + v0 * v1 * v2 * v3 * v47 * v8 - v10 * v42 * v5 * v8 + v1 * v2 * v3 * v42 * v5 * v8
           - 2 * v0 * v2 * v3 * v43 * v5 * v8 - 2 * v1 * v2 * v4 * v43 * v6 * v8
           + v12 * v4 * v44 * v6 * v8 - v10 * v46 * v6 * v8 + v1 * v2 * v3 * v46 * v6 * v8
           - v12 * v3 * v47 * v6 * v8 - v2 * v3 * v41 * v5 * v6 * v8 + v1 * v3 * v43 * v5 * v6 * v8
           - v1 * v2 * v44 * v5 * v6 * v8 - v1 * v14 * v45 * v9 + v17 * v46 * v9
           + v0 * v1 * v39 * v5 * v9 - v14 * v40 * v5 * v9 - v12 * v39 * v7 * v9
           + 2 * v0 * v1 * v40 * v7 * v9 - v14 * v41 * v7 * v9 + v12 * v45 * v8 * v9
           - v0 * v1 * v46 * v8 * v9 - v1 * v40 * v5 * v8 * v9 + v0 * v41 * v5 * v8 * v9);

  gradushift[1][0]
      = v51
        * (v12 * v16 * v36 - v17 * v2 * v3 * v36 - v0 * v1 * v16 * v37 + v10 * v18 * v37
           + 2 * v14 * v15 * v2 * v37 - v1 * v18 * v2 * v3 * v37 - v1 * v11 * v18 * v38
           + v12 * v18 * v3 * v38 + v20 * v3 * v38 - v11 * v14 * v30 * v4 - v12 * v15 * v30 * v4
           - v14 * v15 * v32 * v4 - v11 * v18 * v32 * v4 + 2 * v17 * v2 * v33 * v4
           + 2 * v1 * v18 * v2 * v34 * v4 - v12 * v18 * v35 * v4 - v20 * v35 * v4 + v14 * v38 * v48
           + 2 * v0 * v31 * v4 * v48 - v1 * v16 * v31 * v5 + v0 * v16 * v32 * v5
           + v18 * v2 * v3 * v32 * v5 - v17 * v3 * v33 * v5 - 2 * v14 * v15 * v34 * v5
           - v11 * v18 * v34 * v5 - v1 * v18 * v3 * v34 * v5 + v1 * v18 * v2 * v35 * v5
           + v0 * v33 * v48 * v5 + v2 * v30 * v49 * v5 - 3 * v0 * v2 * v36 * v50 + v14 * v38 * v50
           + 2 * v0 * v31 * v4 * v50 + v2 * v30 * v5 * v50 + v0 * v33 * v5 * v50
           + 2 * v11 * v14 * v36 * v6 - v17 * v3 * v37 * v6 - 2 * v17 * v2 * v38 * v6
           + 2 * v0 * v1 * v2 * v30 * v4 * v6 - 2 * v14 * v2 * v31 * v4 * v6
           - 2 * v1 * v14 * v33 * v4 * v6 + 2 * v17 * v34 * v4 * v6 + v1 * v36 * v48 * v6
           - 3 * v0 * v37 * v48 * v6 + v1 * v36 * v49 * v6 - v0 * v11 * v30 * v5 * v6
           - v0 * v1 * v3 * v30 * v5 * v6 + v17 * v35 * v5 * v6 + v31 * v48 * v5 * v6
           + v31 * v49 * v5 * v6 + v0 * v10 * v30 * v7 - v0 * v1 * v2 * v3 * v30 * v7
           - v0 * v15 * v2 * v32 * v7 - 2 * v11 * v14 * v33 * v7 - v12 * v15 * v33 * v7
           - v17 * v3 * v34 * v7 + v17 * v2 * v35 * v7 + v1 * v33 * v48 * v7 + v0 * v34 * v48 * v7
           + v2 * v31 * v49 * v7 + v1 * v33 * v49 * v7 + v2 * v31 * v50 * v7 + v0 * v34 * v50 * v7
           - v1 * v11 * v30 * v6 * v7 + v12 * v3 * v30 * v6 * v7 - 2 * v0 * v1 * v3 * v31 * v6 * v7
           + 2 * v0 * v1 * v2 * v33 * v6 * v7 - v1 * v14 * v35 * v6 * v7 + v32 * v48 * v6 * v7
           + v32 * v49 * v6 * v7 - v0 * v10 * v36 * v8 + v0 * v1 * v2 * v3 * v36 * v8
           - v12 * v15 * v38 * v8 - 2 * v0 * v1 * v2 * v33 * v4 * v8 - 2 * v14 * v2 * v34 * v4 * v8
           + 2 * v1 * v14 * v35 * v4 * v8 + v2 * v37 * v49 * v8 - 2 * v1 * v38 * v49 * v8
           - 2 * v0 * v2 * v3 * v31 * v5 * v8 - v15 * v2 * v32 * v5 * v8
           + v0 * v1 * v3 * v33 * v5 * v8 - v14 * v2 * v35 * v5 * v8 + v34 * v49 * v5 * v8
           + v2 * v37 * v50 * v8 + v34 * v5 * v50 * v8 - v12 * v3 * v36 * v6 * v8
           + v0 * v1 * v3 * v37 * v6 * v8 + 2 * v0 * v1 * v2 * v38 * v6 * v8
           - 2 * v1 * v2 * v31 * v4 * v6 * v8 + 2 * v0 * v2 * v32 * v4 * v6 * v8
           + 2 * v12 * v33 * v4 * v6 * v8 - 2 * v0 * v1 * v34 * v4 * v6 * v8
           + v1 * v3 * v31 * v5 * v6 * v8 - v0 * v3 * v32 * v5 * v6 * v8
           - 2 * v1 * v2 * v33 * v5 * v6 * v8 + 2 * v0 * v2 * v34 * v5 * v6 * v8
           - v0 * v1 * v35 * v5 * v6 * v8 - v10 * v31 * v7 * v8 + v1 * v2 * v3 * v31 * v7 * v8
           - v0 * v2 * v3 * v32 * v7 * v8 - v12 * v3 * v33 * v7 * v8 + v0 * v1 * v3 * v34 * v7 * v8
           - v0 * v1 * v2 * v35 * v7 * v8 - 2 * v1 * v2 * v34 * v6 * v7 * v8
           + v12 * v35 * v6 * v7 * v8);

  gradushift[1][1]
      = v51
        * (2 * v0 * v11 * v15 * v36 - v1 * v16 * v2 * v36 + v1 * v19 * v37 - 2 * v0 * v16 * v2 * v37
           - v11 * v18 * v3 * v37 - v0 * v1 * v16 * v38 + v10 * v18 * v38 + 2 * v14 * v15 * v2 * v38
           - v1 * v18 * v2 * v3 * v38 - v1 * v16 * v31 * v4 + v0 * v16 * v32 * v4
           + v18 * v2 * v3 * v32 * v4 - v17 * v3 * v33 * v4 - 2 * v14 * v15 * v34 * v4
           - v11 * v18 * v34 * v4 - v1 * v18 * v3 * v34 * v4 + v1 * v18 * v2 * v35 * v4
           + v0 * v3 * v36 * v48 + v15 * v37 * v48 + v0 * v33 * v4 * v48 + v15 * v37 * v49
           + v2 * v30 * v4 * v49 - v11 * v15 * v30 * v5 + 2 * v16 * v2 * v31 * v5
           - 2 * v0 * v15 * v3 * v31 * v5 - v19 * v32 * v5 - 2 * v0 * v15 * v2 * v33 * v5
           + 2 * v0 * v16 * v34 * v5 + 2 * v18 * v2 * v3 * v34 * v5 - v14 * v15 * v35 * v5
           - v11 * v18 * v35 * v5 + v0 * v3 * v36 * v50 + v2 * v30 * v4 * v50 + v0 * v33 * v4 * v50
           - v17 * v3 * v38 * v6 - v0 * v11 * v30 * v4 * v6 - v0 * v1 * v3 * v30 * v4 * v6
           + v17 * v35 * v4 * v6 - 3 * v0 * v38 * v48 * v6 + v31 * v4 * v48 * v6
           - 3 * v2 * v36 * v49 * v6 + v31 * v4 * v49 * v6 + 2 * v0 * v2 * v3 * v30 * v5 * v6
           + 2 * v33 * v48 * v5 * v6 + 2 * v33 * v49 * v5 * v6 - v0 * v11 * v3 * v30 * v7
           - 2 * v11 * v15 * v31 * v7 + v16 * v2 * v32 * v7 - v0 * v15 * v3 * v32 * v7
           - v1 * v16 * v34 * v7 + v3 * v31 * v48 * v7 + v0 * v35 * v48 * v7 + v2 * v33 * v49 * v7
           + v3 * v31 * v50 * v7 + v2 * v33 * v50 * v7 + v0 * v35 * v50 * v7 + v10 * v30 * v6 * v7
           - v1 * v2 * v3 * v30 * v6 * v7 + 2 * v0 * v2 * v3 * v31 * v6 * v7
           - 2 * v0 * v1 * v3 * v33 * v6 * v7 - v14 * v2 * v35 * v6 * v7 + v34 * v48 * v6 * v7
           + v34 * v49 * v6 * v7 - 2 * v0 * v2 * v3 * v31 * v4 * v8 - v15 * v2 * v32 * v4 * v8
           + v0 * v1 * v3 * v33 * v4 * v8 - v14 * v2 * v35 * v4 * v8 + v2 * v38 * v49 * v8
           + v34 * v4 * v49 * v8 + 2 * v15 * v3 * v32 * v5 * v8 - 2 * v0 * v2 * v3 * v33 * v5 * v8
           - 2 * v15 * v2 * v34 * v5 * v8 - 2 * v3 * v37 * v50 * v8 + v2 * v38 * v50 * v8
           + v34 * v4 * v50 * v8 - v10 * v36 * v6 * v8 + v1 * v2 * v3 * v36 * v6 * v8
           + 2 * v0 * v2 * v3 * v37 * v6 * v8 + v0 * v1 * v3 * v38 * v6 * v8
           + v1 * v3 * v31 * v4 * v6 * v8 - v0 * v3 * v32 * v4 * v6 * v8
           - 2 * v1 * v2 * v33 * v4 * v6 * v8 + 2 * v0 * v2 * v34 * v4 * v6 * v8
           - v0 * v1 * v35 * v4 * v6 * v8 - 2 * v2 * v3 * v31 * v5 * v6 * v8
           - 2 * v0 * v3 * v34 * v5 * v6 * v8 + 2 * v0 * v2 * v35 * v5 * v6 * v8
           - v10 * v33 * v7 * v8 + v1 * v2 * v3 * v33 * v7 * v8 - 2 * v0 * v2 * v3 * v34 * v7 * v8
           - v2 * v3 * v32 * v6 * v7 * v8 + v1 * v3 * v34 * v6 * v7 * v8
           - v1 * v2 * v35 * v6 * v7 * v8 + v17 * v36 * v9 + v1 * v18 * v37 * v9
           - v14 * v30 * v5 * v9 - v18 * v32 * v5 * v9 + v0 * v1 * v30 * v7 * v9
           - v14 * v31 * v7 * v9 - v0 * v1 * v36 * v8 * v9 - v14 * v37 * v8 * v9
           + 2 * v0 * v31 * v5 * v8 * v9 - v1 * v31 * v7 * v8 * v9 + v0 * v32 * v7 * v8 * v9);

  gradushift[1][2]
      = v51
        * (-(v12 * v15 * v3 * v36) + 2 * v0 * v11 * v15 * v37 - v1 * v16 * v2 * v37
           + v12 * v16 * v38 - v17 * v2 * v3 * v38 + v0 * v10 * v30 * v4
           - v0 * v1 * v2 * v3 * v30 * v4 - v0 * v15 * v2 * v32 * v4 - 2 * v11 * v14 * v33 * v4
           - v12 * v15 * v33 * v4 - v17 * v3 * v34 * v4 + v17 * v2 * v35 * v4
           - 2 * v1 * v3 * v36 * v48 + v0 * v3 * v37 * v48 + v1 * v33 * v4 * v48
           + v0 * v34 * v4 * v48 + v11 * v36 * v49 + v2 * v31 * v4 * v49 + v1 * v33 * v4 * v49
           - v0 * v11 * v3 * v30 * v5 - 2 * v11 * v15 * v31 * v5 + v16 * v2 * v32 * v5
           - v0 * v15 * v3 * v32 * v5 - v1 * v16 * v34 * v5 + v3 * v31 * v48 * v5
           + v0 * v35 * v48 * v5 + v2 * v33 * v49 * v5 + v11 * v36 * v50 + v0 * v3 * v37 * v50
           - 3 * v0 * v2 * v38 * v50 + v2 * v31 * v4 * v50 + v0 * v34 * v4 * v50
           + v3 * v31 * v5 * v50 + v2 * v33 * v5 * v50 + v0 * v35 * v5 * v50
           - 2 * v0 * v10 * v36 * v6 + 2 * v0 * v1 * v2 * v3 * v36 * v6 + 2 * v11 * v14 * v38 * v6
           - v1 * v11 * v30 * v4 * v6 + v12 * v3 * v30 * v4 * v6 - 2 * v0 * v1 * v3 * v31 * v4 * v6
           + 2 * v0 * v1 * v2 * v33 * v4 * v6 - v1 * v14 * v35 * v4 * v6 + v1 * v38 * v48 * v6
           + v32 * v4 * v48 * v6 - 3 * v2 * v37 * v49 * v6 + v1 * v38 * v49 * v6
           + v32 * v4 * v49 * v6 + v10 * v30 * v5 * v6 - v1 * v2 * v3 * v30 * v5 * v6
           + 2 * v0 * v2 * v3 * v31 * v5 * v6 - 2 * v0 * v1 * v3 * v33 * v5 * v6
           - v14 * v2 * v35 * v5 * v6 + v34 * v48 * v5 * v6 + v34 * v49 * v5 * v6 - v13 * v30 * v7
           + 2 * v1 * v11 * v3 * v30 * v7 - 2 * v0 * v11 * v3 * v31 * v7 - v11 * v15 * v32 * v7
           + 2 * v0 * v10 * v33 * v7 - 2 * v0 * v1 * v2 * v3 * v33 * v7 - v11 * v14 * v35 * v7
           - v12 * v15 * v35 * v7 + 2 * v2 * v34 * v49 * v7 + 2 * v2 * v34 * v50 * v7
           + 2 * v10 * v31 * v6 * v7 - 2 * v1 * v2 * v3 * v31 * v6 * v7
           + 2 * v0 * v2 * v3 * v32 * v6 * v7 - 2 * v1 * v11 * v33 * v6 * v7
           + 2 * v12 * v3 * v33 * v6 * v7 - 2 * v0 * v11 * v34 * v6 * v7
           - 2 * v0 * v1 * v3 * v34 * v6 * v7 + 2 * v0 * v1 * v2 * v35 * v6 * v7 + v13 * v36 * v8
           - v0 * v10 * v38 * v8 + v0 * v1 * v2 * v3 * v38 * v8 - v10 * v31 * v4 * v8
           + v1 * v2 * v3 * v31 * v4 * v8 - v0 * v2 * v3 * v32 * v4 * v8 - v12 * v3 * v33 * v4 * v8
           + v0 * v1 * v3 * v34 * v4 * v8 - v0 * v1 * v2 * v35 * v4 * v8 - v10 * v33 * v5 * v8
           + v1 * v2 * v3 * v33 * v5 * v8 - 2 * v0 * v2 * v3 * v34 * v5 * v8 - v10 * v37 * v6 * v8
           + v1 * v2 * v3 * v37 * v6 * v8 - v12 * v3 * v38 * v6 * v8
           - 2 * v1 * v2 * v34 * v4 * v6 * v8 + v12 * v35 * v4 * v6 * v8
           - v2 * v3 * v32 * v5 * v6 * v8 + v1 * v3 * v34 * v5 * v6 * v8
           - v1 * v2 * v35 * v5 * v6 * v8 - v1 * v14 * v36 * v9 + v17 * v37 * v9
           + v0 * v1 * v30 * v5 * v9 - v14 * v31 * v5 * v9 - v12 * v30 * v7 * v9
           + 2 * v0 * v1 * v31 * v7 * v9 - v14 * v32 * v7 * v9 + v12 * v36 * v8 * v9
           - v0 * v1 * v37 * v8 * v9 - v1 * v31 * v5 * v8 * v9 + v0 * v32 * v5 * v8 * v9);

  gradushift[2][0]
      = v51
        * (v12 * v16 * v27 - v0 * v1 * v16 * v28 + v10 * v18 * v28 + 2 * v14 * v15 * v2 * v28
           - v1 * v11 * v18 * v29 - v17 * v2 * v27 * v3 - v1 * v18 * v2 * v28 * v3
           + v12 * v18 * v29 * v3 + v20 * v29 * v3 - v11 * v14 * v21 * v4 - v12 * v15 * v21 * v4
           - v14 * v15 * v23 * v4 - v11 * v18 * v23 * v4 + 2 * v17 * v2 * v24 * v4
           + 2 * v1 * v18 * v2 * v25 * v4 - v12 * v18 * v26 * v4 - v20 * v26 * v4 + v14 * v29 * v48
           + 2 * v0 * v22 * v4 * v48 - v1 * v16 * v22 * v5 + v0 * v16 * v23 * v5
           - 2 * v14 * v15 * v25 * v5 - v11 * v18 * v25 * v5 + v1 * v18 * v2 * v26 * v5
           + v18 * v2 * v23 * v3 * v5 - v17 * v24 * v3 * v5 - v1 * v18 * v25 * v3 * v5
           + v0 * v24 * v48 * v5 + v2 * v21 * v49 * v5 - 3 * v0 * v2 * v27 * v50 + v14 * v29 * v50
           + 2 * v0 * v22 * v4 * v50 + v2 * v21 * v5 * v50 + v0 * v24 * v5 * v50
           + 2 * v11 * v14 * v27 * v6 - 2 * v17 * v2 * v29 * v6 - v17 * v28 * v3 * v6
           + 2 * v0 * v1 * v2 * v21 * v4 * v6 - 2 * v14 * v2 * v22 * v4 * v6
           - 2 * v1 * v14 * v24 * v4 * v6 + 2 * v17 * v25 * v4 * v6 + v1 * v27 * v48 * v6
           - 3 * v0 * v28 * v48 * v6 + v1 * v27 * v49 * v6 - v0 * v11 * v21 * v5 * v6
           + v17 * v26 * v5 * v6 - v0 * v1 * v21 * v3 * v5 * v6 + v22 * v48 * v5 * v6
           + v22 * v49 * v5 * v6 + v0 * v10 * v21 * v7 - v0 * v15 * v2 * v23 * v7
           - 2 * v11 * v14 * v24 * v7 - v12 * v15 * v24 * v7 + v17 * v2 * v26 * v7
           - v0 * v1 * v2 * v21 * v3 * v7 - v17 * v25 * v3 * v7 + v1 * v24 * v48 * v7
           + v0 * v25 * v48 * v7 + v2 * v22 * v49 * v7 + v1 * v24 * v49 * v7 + v2 * v22 * v50 * v7
           + v0 * v25 * v50 * v7 - v1 * v11 * v21 * v6 * v7 + 2 * v0 * v1 * v2 * v24 * v6 * v7
           - v1 * v14 * v26 * v6 * v7 + v12 * v21 * v3 * v6 * v7 - 2 * v0 * v1 * v22 * v3 * v6 * v7
           + v23 * v48 * v6 * v7 + v23 * v49 * v6 * v7 - v0 * v10 * v27 * v8 - v12 * v15 * v29 * v8
           + v0 * v1 * v2 * v27 * v3 * v8 - 2 * v0 * v1 * v2 * v24 * v4 * v8
           - 2 * v14 * v2 * v25 * v4 * v8 + 2 * v1 * v14 * v26 * v4 * v8 + v2 * v28 * v49 * v8
           - 2 * v1 * v29 * v49 * v8 - v15 * v2 * v23 * v5 * v8 - v14 * v2 * v26 * v5 * v8
           - 2 * v0 * v2 * v22 * v3 * v5 * v8 + v0 * v1 * v24 * v3 * v5 * v8 + v25 * v49 * v5 * v8
           + v2 * v28 * v50 * v8 + v25 * v5 * v50 * v8 + 2 * v0 * v1 * v2 * v29 * v6 * v8
           - v12 * v27 * v3 * v6 * v8 + v0 * v1 * v28 * v3 * v6 * v8
           - 2 * v1 * v2 * v22 * v4 * v6 * v8 + 2 * v0 * v2 * v23 * v4 * v6 * v8
           + 2 * v12 * v24 * v4 * v6 * v8 - 2 * v0 * v1 * v25 * v4 * v6 * v8
           - 2 * v1 * v2 * v24 * v5 * v6 * v8 + 2 * v0 * v2 * v25 * v5 * v6 * v8
           - v0 * v1 * v26 * v5 * v6 * v8 + v1 * v22 * v3 * v5 * v6 * v8
           - v0 * v23 * v3 * v5 * v6 * v8 - v10 * v22 * v7 * v8 - v0 * v1 * v2 * v26 * v7 * v8
           + v1 * v2 * v22 * v3 * v7 * v8 - v0 * v2 * v23 * v3 * v7 * v8 - v12 * v24 * v3 * v7 * v8
           + v0 * v1 * v25 * v3 * v7 * v8 - 2 * v1 * v2 * v25 * v6 * v7 * v8
           + v12 * v26 * v6 * v7 * v8);

  gradushift[2][1]
      = v51
        * (2 * v0 * v11 * v15 * v27 - v1 * v16 * v2 * v27 + v1 * v19 * v28 - 2 * v0 * v16 * v2 * v28
           - v0 * v1 * v16 * v29 + v10 * v18 * v29 + 2 * v14 * v15 * v2 * v29 - v11 * v18 * v28 * v3
           - v1 * v18 * v2 * v29 * v3 - v1 * v16 * v22 * v4 + v0 * v16 * v23 * v4
           - 2 * v14 * v15 * v25 * v4 - v11 * v18 * v25 * v4 + v1 * v18 * v2 * v26 * v4
           + v18 * v2 * v23 * v3 * v4 - v17 * v24 * v3 * v4 - v1 * v18 * v25 * v3 * v4
           + v15 * v28 * v48 + v0 * v27 * v3 * v48 + v0 * v24 * v4 * v48 + v15 * v28 * v49
           + v2 * v21 * v4 * v49 - v11 * v15 * v21 * v5 + 2 * v16 * v2 * v22 * v5 - v19 * v23 * v5
           - 2 * v0 * v15 * v2 * v24 * v5 + 2 * v0 * v16 * v25 * v5 - v14 * v15 * v26 * v5
           - v11 * v18 * v26 * v5 - 2 * v0 * v15 * v22 * v3 * v5 + 2 * v18 * v2 * v25 * v3 * v5
           + v0 * v27 * v3 * v50 + v2 * v21 * v4 * v50 + v0 * v24 * v4 * v50 - v17 * v29 * v3 * v6
           - v0 * v11 * v21 * v4 * v6 + v17 * v26 * v4 * v6 - v0 * v1 * v21 * v3 * v4 * v6
           - 3 * v0 * v29 * v48 * v6 + v22 * v4 * v48 * v6 - 3 * v2 * v27 * v49 * v6
           + v22 * v4 * v49 * v6 + 2 * v0 * v2 * v21 * v3 * v5 * v6 + 2 * v24 * v48 * v5 * v6
           + 2 * v24 * v49 * v5 * v6 - 2 * v11 * v15 * v22 * v7 + v16 * v2 * v23 * v7
           - v1 * v16 * v25 * v7 - v0 * v11 * v21 * v3 * v7 - v0 * v15 * v23 * v3 * v7
           + v0 * v26 * v48 * v7 + v22 * v3 * v48 * v7 + v2 * v24 * v49 * v7 + v2 * v24 * v50 * v7
           + v0 * v26 * v50 * v7 + v22 * v3 * v50 * v7 + v10 * v21 * v6 * v7
           - v14 * v2 * v26 * v6 * v7 - v1 * v2 * v21 * v3 * v6 * v7
           + 2 * v0 * v2 * v22 * v3 * v6 * v7 - 2 * v0 * v1 * v24 * v3 * v6 * v7
           + v25 * v48 * v6 * v7 + v25 * v49 * v6 * v7 - v15 * v2 * v23 * v4 * v8
           - v14 * v2 * v26 * v4 * v8 - 2 * v0 * v2 * v22 * v3 * v4 * v8
           + v0 * v1 * v24 * v3 * v4 * v8 + v2 * v29 * v49 * v8 + v25 * v4 * v49 * v8
           - 2 * v15 * v2 * v25 * v5 * v8 + 2 * v15 * v23 * v3 * v5 * v8
           - 2 * v0 * v2 * v24 * v3 * v5 * v8 + v2 * v29 * v50 * v8 - 2 * v28 * v3 * v50 * v8
           + v25 * v4 * v50 * v8 - v10 * v27 * v6 * v8 + v1 * v2 * v27 * v3 * v6 * v8
           + 2 * v0 * v2 * v28 * v3 * v6 * v8 + v0 * v1 * v29 * v3 * v6 * v8
           - 2 * v1 * v2 * v24 * v4 * v6 * v8 + 2 * v0 * v2 * v25 * v4 * v6 * v8
           - v0 * v1 * v26 * v4 * v6 * v8 + v1 * v22 * v3 * v4 * v6 * v8
           - v0 * v23 * v3 * v4 * v6 * v8 + 2 * v0 * v2 * v26 * v5 * v6 * v8
           - 2 * v2 * v22 * v3 * v5 * v6 * v8 - 2 * v0 * v25 * v3 * v5 * v6 * v8
           - v10 * v24 * v7 * v8 + v1 * v2 * v24 * v3 * v7 * v8 - 2 * v0 * v2 * v25 * v3 * v7 * v8
           - v1 * v2 * v26 * v6 * v7 * v8 - v2 * v23 * v3 * v6 * v7 * v8
           + v1 * v25 * v3 * v6 * v7 * v8 + v17 * v27 * v9 + v1 * v18 * v28 * v9
           - v14 * v21 * v5 * v9 - v18 * v23 * v5 * v9 + v0 * v1 * v21 * v7 * v9
           - v14 * v22 * v7 * v9 - v0 * v1 * v27 * v8 * v9 - v14 * v28 * v8 * v9
           + 2 * v0 * v22 * v5 * v8 * v9 - v1 * v22 * v7 * v8 * v9 + v0 * v23 * v7 * v8 * v9);

  gradushift[2][2]
      = v51
        * (2 * v0 * v11 * v15 * v28 - v1 * v16 * v2 * v28 + v12 * v16 * v29 - v12 * v15 * v27 * v3
           - v17 * v2 * v29 * v3 + v0 * v10 * v21 * v4 - v0 * v15 * v2 * v23 * v4
           - 2 * v11 * v14 * v24 * v4 - v12 * v15 * v24 * v4 + v17 * v2 * v26 * v4
           - v0 * v1 * v2 * v21 * v3 * v4 - v17 * v25 * v3 * v4 - 2 * v1 * v27 * v3 * v48
           + v0 * v28 * v3 * v48 + v1 * v24 * v4 * v48 + v0 * v25 * v4 * v48 + v11 * v27 * v49
           + v2 * v22 * v4 * v49 + v1 * v24 * v4 * v49 - 2 * v11 * v15 * v22 * v5
           + v16 * v2 * v23 * v5 - v1 * v16 * v25 * v5 - v0 * v11 * v21 * v3 * v5
           - v0 * v15 * v23 * v3 * v5 + v0 * v26 * v48 * v5 + v22 * v3 * v48 * v5
           + v2 * v24 * v49 * v5 + v11 * v27 * v50 - 3 * v0 * v2 * v29 * v50 + v0 * v28 * v3 * v50
           + v2 * v22 * v4 * v50 + v0 * v25 * v4 * v50 + v2 * v24 * v5 * v50 + v0 * v26 * v5 * v50
           + v22 * v3 * v5 * v50 - 2 * v0 * v10 * v27 * v6 + 2 * v11 * v14 * v29 * v6
           + 2 * v0 * v1 * v2 * v27 * v3 * v6 - v1 * v11 * v21 * v4 * v6
           + 2 * v0 * v1 * v2 * v24 * v4 * v6 - v1 * v14 * v26 * v4 * v6 + v12 * v21 * v3 * v4 * v6
           - 2 * v0 * v1 * v22 * v3 * v4 * v6 + v1 * v29 * v48 * v6 + v23 * v4 * v48 * v6
           - 3 * v2 * v28 * v49 * v6 + v1 * v29 * v49 * v6 + v23 * v4 * v49 * v6
           + v10 * v21 * v5 * v6 - v14 * v2 * v26 * v5 * v6 - v1 * v2 * v21 * v3 * v5 * v6
           + 2 * v0 * v2 * v22 * v3 * v5 * v6 - 2 * v0 * v1 * v24 * v3 * v5 * v6
           + v25 * v48 * v5 * v6 + v25 * v49 * v5 * v6 - v13 * v21 * v7 - v11 * v15 * v23 * v7
           + 2 * v0 * v10 * v24 * v7 - v11 * v14 * v26 * v7 - v12 * v15 * v26 * v7
           + 2 * v1 * v11 * v21 * v3 * v7 - 2 * v0 * v11 * v22 * v3 * v7
           - 2 * v0 * v1 * v2 * v24 * v3 * v7 + 2 * v2 * v25 * v49 * v7 + 2 * v2 * v25 * v50 * v7
           + 2 * v10 * v22 * v6 * v7 - 2 * v1 * v11 * v24 * v6 * v7 - 2 * v0 * v11 * v25 * v6 * v7
           + 2 * v0 * v1 * v2 * v26 * v6 * v7 - 2 * v1 * v2 * v22 * v3 * v6 * v7
           + 2 * v0 * v2 * v23 * v3 * v6 * v7 + 2 * v12 * v24 * v3 * v6 * v7
           - 2 * v0 * v1 * v25 * v3 * v6 * v7 + v13 * v27 * v8 - v0 * v10 * v29 * v8
           + v0 * v1 * v2 * v29 * v3 * v8 - v10 * v22 * v4 * v8 - v0 * v1 * v2 * v26 * v4 * v8
           + v1 * v2 * v22 * v3 * v4 * v8 - v0 * v2 * v23 * v3 * v4 * v8 - v12 * v24 * v3 * v4 * v8
           + v0 * v1 * v25 * v3 * v4 * v8 - v10 * v24 * v5 * v8 + v1 * v2 * v24 * v3 * v5 * v8
           - 2 * v0 * v2 * v25 * v3 * v5 * v8 - v10 * v28 * v6 * v8 + v1 * v2 * v28 * v3 * v6 * v8
           - v12 * v29 * v3 * v6 * v8 - 2 * v1 * v2 * v25 * v4 * v6 * v8 + v12 * v26 * v4 * v6 * v8
           - v1 * v2 * v26 * v5 * v6 * v8 - v2 * v23 * v3 * v5 * v6 * v8
           + v1 * v25 * v3 * v5 * v6 * v8 - v1 * v14 * v27 * v9 + v17 * v28 * v9
           + v0 * v1 * v21 * v5 * v9 - v14 * v22 * v5 * v9 - v12 * v21 * v7 * v9
           + 2 * v0 * v1 * v22 * v7 * v9 - v14 * v23 * v7 * v9 + v12 * v27 * v8 * v9
           - v0 * v1 * v28 * v8 * v9 - v1 * v22 * v5 * v8 * v9 + v0 * v23 * v5 * v8 * v9);

  return gradushift;
}