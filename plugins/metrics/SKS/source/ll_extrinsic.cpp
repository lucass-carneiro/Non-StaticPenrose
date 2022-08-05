#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::ll_extrinsic(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {

  using std::pow;

  metric_server::spatial_matrix llK{};

  const double v0 = llgSKS_33(t, x, y, z);
  const double v1 = llgSKS_22(t, x, y, z);
  const double v2 = llgSKS_11(t, x, y, z);
  const double v3 = llgSKS_12(t, x, y, z);
  const double v4 = llgSKS_23(t, x, y, z);
  const double v5 = llgSKS_13(t, x, y, z);
  const double v6 = llgSKS_01(t, x, y, z);
  const double v7 = llgSKS_02(t, x, y, z);
  const double v8 = llgSKS_03(t, x, y, z);
  const double v9 = llgSKS_00(t, x, y, z);
  const double v10 = pow(v3, 2);
  const double v11 = pow(v4, 2);
  const double v12 = pow(v5, 2);
  const double v13 = pow(v6, 2);
  const double v14 = pow(v7, 2);
  const double v15 = pow(v8, 2);
  const double v16 = dllgSKS_33_dt(t, x, y, z);
  const double v17 = dllgSKS_33_dx(t, x, y, z);
  const double v18 = dllgSKS_33_dy(t, x, y, z);
  const double v19 = dllgSKS_33_dz(t, x, y, z);
  const double v20 = dllgSKS_23_dz(t, x, y, z);
  const double v21 = dllgSKS_13_dz(t, x, y, z);
  const double v22 = dllgSKS_03_dz(t, x, y, z);
  const double v23 = dllgSKS_23_dt(t, x, y, z);
  const double v24 = dllgSKS_23_dx(t, x, y, z);
  const double v25 = dllgSKS_13_dy(t, x, y, z);
  const double v26 = dllgSKS_03_dy(t, x, y, z);
  const double v27 = dllgSKS_22_dz(t, x, y, z);
  const double v28 = dllgSKS_12_dz(t, x, y, z);
  const double v29 = dllgSKS_02_dz(t, x, y, z);
  const double v30 = dllgSKS_13_dt(t, x, y, z);
  const double v31 = dllgSKS_03_dx(t, x, y, z);
  const double v32 = dllgSKS_11_dz(t, x, y, z);
  const double v33 = dllgSKS_01_dz(t, x, y, z);
  const double v34 = dllgSKS_22_dt(t, x, y, z);
  const double v35 = dllgSKS_22_dx(t, x, y, z);
  const double v36 = dllgSKS_23_dy(t, x, y, z);
  const double v37 = dllgSKS_22_dy(t, x, y, z);
  const double v38 = dllgSKS_12_dy(t, x, y, z);
  const double v39 = dllgSKS_02_dy(t, x, y, z);
  const double v40 = dllgSKS_12_dt(t, x, y, z);
  const double v41 = dllgSKS_02_dx(t, x, y, z);
  const double v42 = dllgSKS_11_dy(t, x, y, z);
  const double v43 = dllgSKS_01_dy(t, x, y, z);
  const double v44 = dllgSKS_11_dt(t, x, y, z);
  const double v45 = dllgSKS_13_dx(t, x, y, z);
  const double v46 = dllgSKS_12_dx(t, x, y, z);
  const double v47 = dllgSKS_11_dx(t, x, y, z);
  const double v48 = dllgSKS_01_dx(t, x, y, z);
  const double v49 = v0 * v10;
  const double v50 = v11 * v2;
  const double v51 = v1 * v12;
  const double v52 = -(v0 * v1 * v2);
  const double v53 = -2 * v3 * v4 * v5;
  const double v54 = 1 / (v49 + v50 + v51 + v52 + v53);
  const double v55 = 1
                     / sqrt(v54
                            * (-(v0 * v1 * v13) + v11 * v13 + v12 * v14 + v10 * v15 - v0 * v14 * v2
                               - v1 * v15 * v2 + 2 * v0 * v3 * v6 * v7 - 2 * v4 * v5 * v6 * v7
                               - 2 * v3 * v4 * v6 * v8 + 2 * v1 * v5 * v6 * v8
                               + 2 * v2 * v4 * v7 * v8 - 2 * v3 * v5 * v7 * v8 + v0 * v1 * v2 * v9
                               - v49 * v9 + 2 * v3 * v4 * v5 * v9 - v50 * v9 - v51 * v9));
  const double v56
      = (v54 * v55
         * (v0 * v1 * v2 * v23 - v23 * v49 + v26 * v49 + v29 * v49 + 2 * v23 * v3 * v4 * v5
            - v23 * v50 + v26 * v50 + v29 * v50 - v23 * v51 + v26 * v51 + v29 * v51 + v26 * v52
            + v29 * v52 + v26 * v53 + v29 * v53 - v0 * v1 * v24 * v6 + v11 * v24 * v6
            + v0 * v1 * v25 * v6 - v11 * v25 * v6 + v0 * v1 * v28 * v6 - v11 * v28 * v6
            - v0 * v27 * v3 * v6 + v18 * v3 * v4 * v6 - v1 * v18 * v5 * v6 + v27 * v4 * v5 * v6
            - v12 * v27 * v7 + v0 * v2 * v27 * v7 + v0 * v24 * v3 * v7 - v0 * v25 * v3 * v7
            - v0 * v28 * v3 * v7 - v18 * v2 * v4 * v7 + v18 * v3 * v5 * v7 - v24 * v4 * v5 * v7
            + v25 * v4 * v5 * v7 + v28 * v4 * v5 * v7 - v10 * v18 * v8 + v1 * v18 * v2 * v8
            - v2 * v27 * v4 * v8 - v24 * v3 * v4 * v8 + v25 * v3 * v4 * v8 + v28 * v3 * v4 * v8
            + v1 * v24 * v5 * v8 - v1 * v25 * v5 * v8 - v1 * v28 * v5 * v8 + v27 * v3 * v5 * v8))
        / 2.;
  const double v57
      = (v54 * v55
         * (v0 * v1 * v2 * v30 - v30 * v49 + v31 * v49 + v33 * v49 + 2 * v3 * v30 * v4 * v5
            - v30 * v50 + v31 * v50 + v33 * v50 - v30 * v51 + v31 * v51 + v33 * v51 + v31 * v52
            + v33 * v52 + v31 * v53 + v33 * v53 - v0 * v24 * v3 * v6 + v0 * v25 * v3 * v6
            - v0 * v28 * v3 * v6 + v0 * v1 * v32 * v6 - v11 * v32 * v6 + v17 * v3 * v4 * v6
            - v1 * v17 * v5 * v6 + v24 * v4 * v5 * v6 - v25 * v4 * v5 * v6 + v28 * v4 * v5 * v6
            - v12 * v24 * v7 + v0 * v2 * v24 * v7 + v12 * v25 * v7 - v0 * v2 * v25 * v7
            - v12 * v28 * v7 + v0 * v2 * v28 * v7 - v0 * v3 * v32 * v7 - v17 * v2 * v4 * v7
            + v17 * v3 * v5 * v7 + v32 * v4 * v5 * v7 - v10 * v17 * v8 + v1 * v17 * v2 * v8
            - v2 * v24 * v4 * v8 + v2 * v25 * v4 * v8 - v2 * v28 * v4 * v8 + v3 * v32 * v4 * v8
            + v24 * v3 * v5 * v8 - v25 * v3 * v5 * v8 + v28 * v3 * v5 * v8 - v1 * v32 * v5 * v8))
        / 2.;
  const double v58
      = (v54 * v55
         * (v0 * v1 * v2 * v40 - v40 * v49 + v41 * v49 + v43 * v49 + 2 * v3 * v4 * v40 * v5
            - v40 * v50 + v41 * v50 + v43 * v50 - v40 * v51 + v41 * v51 + v43 * v51 + v41 * v52
            + v43 * v52 + v41 * v53 + v43 * v53 - v0 * v3 * v35 * v6 + v24 * v3 * v4 * v6
            + v25 * v3 * v4 * v6 - v28 * v3 * v4 * v6 + v0 * v1 * v42 * v6 - v11 * v42 * v6
            - v1 * v24 * v5 * v6 - v1 * v25 * v5 * v6 + v1 * v28 * v5 * v6 + v35 * v4 * v5 * v6
            - v12 * v35 * v7 + v0 * v2 * v35 * v7 - v2 * v24 * v4 * v7 - v2 * v25 * v4 * v7
            + v2 * v28 * v4 * v7 - v0 * v3 * v42 * v7 + v24 * v3 * v5 * v7 + v25 * v3 * v5 * v7
            - v28 * v3 * v5 * v7 + v4 * v42 * v5 * v7 - v10 * v24 * v8 + v1 * v2 * v24 * v8
            - v10 * v25 * v8 + v1 * v2 * v25 * v8 + v10 * v28 * v8 - v1 * v2 * v28 * v8
            - v2 * v35 * v4 * v8 + v3 * v4 * v42 * v8 + v3 * v35 * v5 * v8 - v1 * v42 * v5 * v8))
        / 2.;

  llK[0][0]
      = (v54 * v55
         * (v0 * v1 * v2 * v44 - 2 * v0 * v1 * v2 * v48 - v44 * v49 + 2 * v48 * v49
            + 2 * v3 * v4 * v44 * v5 - 4 * v3 * v4 * v48 * v5 - v44 * v50 + 2 * v48 * v50
            - v44 * v51 + 2 * v48 * v51 - v3 * v32 * v4 * v6 + v0 * v3 * v42 * v6
            + 2 * v3 * v4 * v45 * v6 - 2 * v0 * v3 * v46 * v6 + v0 * v1 * v47 * v6 - v11 * v47 * v6
            + v1 * v32 * v5 * v6 - v4 * v42 * v5 * v6 - 2 * v1 * v45 * v5 * v6
            + 2 * v4 * v46 * v5 * v6 + v2 * v32 * v4 * v7 + v12 * v42 * v7 - v0 * v2 * v42 * v7
            - 2 * v2 * v4 * v45 * v7 - 2 * v12 * v46 * v7 + 2 * v0 * v2 * v46 * v7
            - v0 * v3 * v47 * v7 - v3 * v32 * v5 * v7 + 2 * v3 * v45 * v5 * v7 + v4 * v47 * v5 * v7
            + v10 * v32 * v8 - v1 * v2 * v32 * v8 + v2 * v4 * v42 * v8 - 2 * v10 * v45 * v8
            + 2 * v1 * v2 * v45 * v8 - 2 * v2 * v4 * v46 * v8 + v3 * v4 * v47 * v8
            - v3 * v42 * v5 * v8 + 2 * v3 * v46 * v5 * v8 - v1 * v47 * v5 * v8))
        / 2.;
  llK[0][1] = v58;
  llK[0][2] = v57;
  llK[1][0] = v58;
  llK[1][1]
      = (v54 * v55
         * (v0 * v1 * v2 * v34 - 2 * v0 * v1 * v2 * v39 - v34 * v49 + 2 * v39 * v49
            + 2 * v3 * v34 * v4 * v5 - 4 * v3 * v39 * v4 * v5 - v34 * v50 + 2 * v39 * v50
            - v34 * v51 + 2 * v39 * v51 - v0 * v1 * v35 * v6 + v11 * v35 * v6 - v0 * v3 * v37 * v6
            + 2 * v0 * v1 * v38 * v6 - 2 * v11 * v38 * v6 - v27 * v3 * v4 * v6
            + 2 * v3 * v36 * v4 * v6 + v1 * v27 * v5 * v6 - 2 * v1 * v36 * v5 * v6
            + v37 * v4 * v5 * v6 + v0 * v3 * v35 * v7 - v12 * v37 * v7 + v0 * v2 * v37 * v7
            - 2 * v0 * v3 * v38 * v7 + v2 * v27 * v4 * v7 - 2 * v2 * v36 * v4 * v7
            - v27 * v3 * v5 * v7 + 2 * v3 * v36 * v5 * v7 - v35 * v4 * v5 * v7
            + 2 * v38 * v4 * v5 * v7 + v10 * v27 * v8 - v1 * v2 * v27 * v8 - 2 * v10 * v36 * v8
            + 2 * v1 * v2 * v36 * v8 - v3 * v35 * v4 * v8 - v2 * v37 * v4 * v8
            + 2 * v3 * v38 * v4 * v8 + v1 * v35 * v5 * v8 + v3 * v37 * v5 * v8
            - 2 * v1 * v38 * v5 * v8))
        / 2.;
  llK[1][2] = v56;
  llK[2][0] = v57;
  llK[2][1] = v56;
  llK[2][2]
      = (v54 * v55
         * (v0 * v1 * v16 * v2 - 2 * v0 * v1 * v2 * v22 - v16 * v49 + 2 * v22 * v49
            + 2 * v16 * v3 * v4 * v5 - 4 * v22 * v3 * v4 * v5 - v16 * v50 + 2 * v22 * v50
            - v16 * v51 + 2 * v22 * v51 - v0 * v1 * v17 * v6 + v11 * v17 * v6
            + 2 * v0 * v1 * v21 * v6 - 2 * v11 * v21 * v6 + v0 * v18 * v3 * v6
            - 2 * v0 * v20 * v3 * v6 + v19 * v3 * v4 * v6 - v1 * v19 * v5 * v6 - v18 * v4 * v5 * v6
            + 2 * v20 * v4 * v5 * v6 + v12 * v18 * v7 - v0 * v18 * v2 * v7 - 2 * v12 * v20 * v7
            + 2 * v0 * v2 * v20 * v7 + v0 * v17 * v3 * v7 - 2 * v0 * v21 * v3 * v7
            - v19 * v2 * v4 * v7 + v19 * v3 * v5 * v7 - v17 * v4 * v5 * v7 + 2 * v21 * v4 * v5 * v7
            - v10 * v19 * v8 + v1 * v19 * v2 * v8 + v18 * v2 * v4 * v8 - 2 * v2 * v20 * v4 * v8
            - v17 * v3 * v4 * v8 + 2 * v21 * v3 * v4 * v8 + v1 * v17 * v5 * v8
            - 2 * v1 * v21 * v5 * v8 - v18 * v3 * v5 * v8 + 2 * v20 * v3 * v5 * v8))
        / 2.;

  return llK;
}