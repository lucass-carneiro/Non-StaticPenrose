#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::spatial_christoffel(double t, double x, double y, double z)
    -> metric_server::chirstofell_t {

  using std::pow;

  metric_server::chirstofell_t Gamma{};

  const double v0 = llgSKS_23(t, x, y, z);
  const double v1 = llgSKS_12(t, x, y, z);
  const double v2 = llgSKS_22(t, x, y, z);
  const double v3 = llgSKS_13(t, x, y, z);
  const double v4 = llgSKS_11(t, x, y, z);
  const double v5 = llgSKS_33(t, x, y, z);
  const double v6 = pow(v1, 2);
  const double v7 = pow(v0, 2);
  const double v8 = pow(v3, 2);
  const double v9 = dllgSKS_33_dx(t, x, y, z);
  const double v10 = dllgSKS_33_dy(t, x, y, z);
  const double v11 = dllgSKS_33_dz(t, x, y, z);
  const double v12 = dllgSKS_23_dz(t, x, y, z);
  const double v13 = dllgSKS_13_dz(t, x, y, z);
  const double v14 = dllgSKS_23_dx(t, x, y, z);
  const double v15 = dllgSKS_13_dy(t, x, y, z);
  const double v16 = dllgSKS_22_dz(t, x, y, z);
  const double v17 = dllgSKS_12_dz(t, x, y, z);
  const double v18 = dllgSKS_11_dz(t, x, y, z);
  const double v19 = dllgSKS_22_dx(t, x, y, z);
  const double v20 = dllgSKS_23_dy(t, x, y, z);
  const double v21 = dllgSKS_22_dy(t, x, y, z);
  const double v22 = dllgSKS_12_dy(t, x, y, z);
  const double v23 = dllgSKS_11_dy(t, x, y, z);
  const double v24 = dllgSKS_13_dx(t, x, y, z);
  const double v25 = dllgSKS_12_dx(t, x, y, z);
  const double v26 = dllgSKS_11_dx(t, x, y, z);
  const double v27 = v0 * v1 * v14;
  const double v28 = v15 * v2 * v3;
  const double v29 = v0 * v16 * v4;
  const double v30 = v17 * v2 * v3;
  const double v31 = v0 * v14 * v4;
  const double v32 = v1 * v15 * v3;
  const double v33 = v0 * v17 * v4;
  const double v34 = v18 * v2 * v3;
  const double v35 = v0 * v14 * v3;
  const double v36 = v0 * v10 * v4;
  const double v37 = v1 * v15 * v5;
  const double v38 = v1 * v17 * v5;
  const double v39 = v0 * v15 * v4;
  const double v40 = v1 * v23 * v5;
  const double v41 = v1 * v17 * v3;
  const double v42 = v2 * v3 * v9;
  const double v43 = v1 * v14 * v5;
  const double v44 = v0 * v15 * v3;
  const double v45 = v14 * v2 * v3;
  const double v46 = v1 * v19 * v5;
  const double v47 = v0 * v1 * v17;
  const double v48 = -(v0 * v1 * v15);
  const double v49 = -(v1 * v14 * v3);
  const double v50 = -(v0 * v17 * v3);
  const double v51 = 1 / (-2 * v0 * v1 * v3 - v2 * v4 * v5 + v5 * v6 + v4 * v7 + v2 * v8);
  const double v52
      = (v51
         * (-(v0 * v1 * v23) - v1 * v19 * v3 + v2 * v23 * v3 + v0 * v19 * v4 - v14 * v2 * v4
            - v15 * v2 * v4 + v17 * v2 * v4 + v14 * v6 + v15 * v6 - v17 * v6))
        / 2.;
  const double v53
      = (v51
         * (-(v0 * v18 * v3) + v1 * v18 * v5 - v14 * v4 * v5 + v15 * v4 * v5 - v17 * v4 * v5
            + v14 * v8 - v15 * v8 + v17 * v8 - v1 * v3 * v9 + v0 * v4 * v9))
        / 2.;
  const double v54
      = (v51
         * (-(v0 * v1 * v10) - v0 * v16 * v3 + v10 * v2 * v3 + v1 * v16 * v5 + v14 * v2 * v5
            - v15 * v2 * v5 - v17 * v2 * v5 - v14 * v7 + v15 * v7 + v17 * v7))
        / 2.;
  const double v55
      = (v51 * (v27 + v28 + v29 - v1 * v16 * v3 + v30 - v10 * v2 * v4 - v45 - v47 + v48 + v10 * v6))
        / 2.;
  const double v56
      = (v51
         * (-(v0 * v1 * v18) + v31 + v32 + v33 + v34 - v39 - v41 + v49 - v2 * v4 * v9 + v6 * v9))
        / 2.;
  const double v57
      = (v51
         * (-(v1 * v10 * v3) + v35 + v36 + v37 + v38 - v43 - v44 - v16 * v4 * v5 + v50 + v16 * v8))
        / 2.;
  const double v58
      = (v51
         * (-(v0 * v23 * v3) + v31 - v32 - v33 + v39 + v40 + v41 + v49 - v19 * v4 * v5 + v19 * v8))
        / 2.;
  const double v59
      = (v51 * (-v35 - v37 + v38 + v42 + v43 + v44 - v18 * v2 * v5 + v50 + v18 * v7 - v0 * v1 * v9))
        / 2.;
  const double v60
      = (v51
         * (-v27 + v28 - v0 * v19 * v3 - v30 + v45 + v46 + v47 + v48 - v2 * v23 * v5 + v23 * v7))
        / 2.;

  Gamma[0][0][0]
      = (v51
         * (v0 * v1 * v18 - 2 * v0 * v1 * v24 + v0 * v23 * v3 + 2 * v2 * v24 * v3
            - 2 * v0 * v25 * v3 - v34 - v40 + 2 * v1 * v25 * v5 - v2 * v26 * v5 + v26 * v7))
        / 2.;

  Gamma[0][0][1] = v60;

  Gamma[0][0][2] = v59;

  Gamma[0][1][0] = v60;

  Gamma[0][1][1]
      = (v51
         * (v0 * v1 * v16 - 2 * v0 * v1 * v20 - v16 * v2 * v3 + 2 * v2 * v20 * v3 - v0 * v21 * v3
            + v19 * v2 * v5 + v1 * v21 * v5 - 2 * v2 * v22 * v5 - v19 * v7 + 2 * v22 * v7))
        / 2.;

  Gamma[0][1][2] = v54;

  Gamma[0][2][0] = v59;

  Gamma[0][2][1] = v54;

  Gamma[0][2][2]
      = (v51
         * (-(v0 * v1 * v11) + v0 * v10 * v3 - 2 * v0 * v12 * v3 + v11 * v2 * v3 - v1 * v10 * v5
            + 2 * v1 * v12 * v5 - 2 * v13 * v2 * v5 + 2 * v13 * v7 + v2 * v5 * v9 - v7 * v9))
        / 2.;

  Gamma[1][0][0]
      = (v51
         * (v1 * v18 * v3 - 2 * v1 * v24 * v3 - v0 * v26 * v3 - v0 * v18 * v4 + 2 * v0 * v24 * v4
            + v1 * v26 * v5 + v23 * v4 * v5 - 2 * v25 * v4 * v5 - v23 * v8 + 2 * v25 * v8))
        / 2.;

  Gamma[1][0][1] = v58;

  Gamma[1][0][2] = v53;

  Gamma[1][1][0] = v58;

  Gamma[1][1][1] = (v51
                    * (-v29 + v1 * v16 * v3 + v0 * v19 * v3 - 2 * v1 * v20 * v3 - 2 * v0 * v22 * v3
                       + 2 * v0 * v20 * v4 - v46 + 2 * v1 * v22 * v5 - v21 * v4 * v5 + v21 * v8))
                   / 2.;

  Gamma[1][1][2] = v57;

  Gamma[1][2][0] = v53;

  Gamma[1][2][1] = v57;

  Gamma[1][2][2]
      = (v51
         * (-(v1 * v11 * v3) - 2 * v0 * v13 * v3 + v0 * v11 * v4 + 2 * v1 * v13 * v5 + v10 * v4 * v5
            - 2 * v12 * v4 * v5 - v10 * v8 + 2 * v12 * v8 + v0 * v3 * v9 - v1 * v5 * v9))
        / 2.;

  Gamma[2][0][0]
      = (v51
         * (-(v0 * v1 * v26) + v1 * v23 * v3 - 2 * v1 * v25 * v3 + v2 * v26 * v3 + v18 * v2 * v4
            - v0 * v23 * v4 - 2 * v2 * v24 * v4 + 2 * v0 * v25 * v4 - v18 * v6 + 2 * v24 * v6))
        / 2.;

  Gamma[2][0][1] = v52;

  Gamma[2][0][2] = v56;

  Gamma[2][1][0] = v52;

  Gamma[2][1][1]
      = (v51
         * (v0 * v1 * v19 - 2 * v0 * v1 * v22 - v19 * v2 * v3 - v1 * v21 * v3 + 2 * v2 * v22 * v3
            + v16 * v2 * v4 - 2 * v2 * v20 * v4 + v0 * v21 * v4 - v16 * v6 + 2 * v20 * v6))
        / 2.;

  Gamma[2][1][2] = v55;

  Gamma[2][2][0] = v56;

  Gamma[2][2][1] = v55;

  Gamma[2][2][2] = (v51
                    * (-2 * v0 * v1 * v13 + v1 * v10 * v3 - 2 * v1 * v12 * v3 + 2 * v13 * v2 * v3
                       - v36 + 2 * v0 * v12 * v4 - v11 * v2 * v4 - v42 + v11 * v6 + v0 * v1 * v9))
                   / 2.;

  return Gamma;
}