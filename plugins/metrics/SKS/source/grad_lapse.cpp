#include "SKS.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::grad_lapse(double t, double x, double y, double z)
    -> metric_server::spatial_vector {

  using std::pow;

  metric_server::spatial_vector gradalpha{};

  const double v0{llgSKS_23(t, x, y, z)};
  const double v1{llgSKS_12(t, x, y, z)};
  const double v2{llgSKS_01(t, x, y, z)};
  const double v3{llgSKS_11(t, x, y, z)};
  const double v4{llgSKS_02(t, x, y, z)};
  const double v5{llgSKS_22(t, x, y, z)};
  const double v6{llgSKS_13(t, x, y, z)};
  const double v7{llgSKS_03(t, x, y, z)};
  const double v8{llgSKS_33(t, x, y, z)};
  const double v9{llgSKS_00(t, x, y, z)};
  const double v10{pow(v0, 2)};
  const double v11{pow(v1, 2)};
  const double v12{pow(v2, 2)};
  const double v13{pow(v3, 2)};
  const double v14{pow(v4, 2)};
  const double v15{pow(v1, 3)};
  const double v16{pow(v5, 2)};
  const double v17{pow(v6, 2)};
  const double v18{pow(v7, 2)};
  const double v19{pow(v1, 4)};
  const double v20{pow(v6, 3)};
  const double v21{pow(v8, 2)};
  const double v22{pow(v6, 4)};
  const double v23{pow(v0, 3)};
  const double v24{pow(v0, 4)};
  const double v25{dllgSKS_33_dz(t, x, y, z)};
  const double v26{dllgSKS_23_dz(t, x, y, z)};
  const double v27{dllgSKS_22_dz(t, x, y, z)};
  const double v28{dllgSKS_13_dz(t, x, y, z)};
  const double v29{dllgSKS_12_dz(t, x, y, z)};
  const double v30{dllgSKS_11_dz(t, x, y, z)};
  const double v31{dllgSKS_03_dz(t, x, y, z)};
  const double v32{dllgSKS_02_dz(t, x, y, z)};
  const double v33{dllgSKS_01_dz(t, x, y, z)};
  const double v34{dllgSKS_00_dz(t, x, y, z)};
  const double v35{dllgSKS_33_dy(t, x, y, z)};
  const double v36{dllgSKS_23_dy(t, x, y, z)};
  const double v37{dllgSKS_22_dy(t, x, y, z)};
  const double v38{dllgSKS_13_dy(t, x, y, z)};
  const double v39{dllgSKS_12_dy(t, x, y, z)};
  const double v40{dllgSKS_11_dy(t, x, y, z)};
  const double v41{dllgSKS_03_dy(t, x, y, z)};
  const double v42{dllgSKS_02_dy(t, x, y, z)};
  const double v43{dllgSKS_01_dy(t, x, y, z)};
  const double v44{dllgSKS_00_dy(t, x, y, z)};
  const double v45{dllgSKS_33_dx(t, x, y, z)};
  const double v46{dllgSKS_23_dx(t, x, y, z)};
  const double v47{dllgSKS_22_dx(t, x, y, z)};
  const double v48{dllgSKS_13_dx(t, x, y, z)};
  const double v49{dllgSKS_12_dx(t, x, y, z)};
  const double v50{dllgSKS_11_dx(t, x, y, z)};
  const double v51{dllgSKS_03_dx(t, x, y, z)};
  const double v52{dllgSKS_02_dx(t, x, y, z)};
  const double v53{dllgSKS_01_dx(t, x, y, z)};
  const double v54{dllgSKS_00_dx(t, x, y, z)};
  const double v55{v11 * v8};
  const double v56{v10 * v3};
  const double v57{v17 * v5};
  const double v58{v10 * v12};
  const double v59{v14 * v17};
  const double v60{v11 * v18};
  const double v61{-2 * v0 * v2 * v4 * v6};
  const double v62{-2 * v0 * v1 * v2 * v7};
  const double v63{-2 * v1 * v4 * v6 * v7};
  const double v64{v55 + v56 + v57 - 2 * v0 * v1 * v6 - v3 * v5 * v8};
  const double v65{pow(v64, -2)};
  const double v66{
      1
      / sqrt((-(v18 * v3 * v5) + v58 + v59 + v60 + v61 + v62 + v63 + 2 * v0 * v3 * v4 * v7
              + 2 * v2 * v5 * v6 * v7 - v14 * v3 * v8 + 2 * v1 * v2 * v4 * v8 - v12 * v5 * v8
              - v55 * v9 - v56 * v9 - v57 * v9 + 2 * v0 * v1 * v6 * v9 + v3 * v5 * v8 * v9)
             / v64)};

  gradalpha[0]
      = (v65 * v66
         * (-(v10 * v13 * v14 * v45) - v12 * v16 * v17 * v45 - v13 * v16 * v18 * v45
            - v18 * v19 * v45 + 2 * v1 * v14 * v20 * v46 - v10 * v13 * v18 * v47
            - v11 * v12 * v21 * v47 - v13 * v14 * v21 * v47 - v14 * v22 * v47
            + 2 * v0 * v2 * v20 * v4 * v47 + 2 * v1 * v2 * v21 * v3 * v4 * v47
            + 2 * v0 * v15 * v18 * v48 + 2 * v1 * v12 * v23 * v48 - 2 * v2 * v23 * v3 * v4 * v48
            + 2 * v0 * v14 * v20 * v49 + 2 * v1 * v14 * v21 * v3 * v49
            - 4 * v10 * v17 * v2 * v4 * v49 - 2 * v11 * v2 * v21 * v4 * v49
            + 2 * v0 * v13 * v18 * v46 * v5 - 2 * v2 * v20 * v4 * v46 * v5
            - 2 * v0 * v1 * v18 * v3 * v48 * v5 + 2 * v1 * v12 * v21 * v49 * v5
            - 2 * v2 * v21 * v3 * v4 * v49 * v5 - v16 * v17 * v18 * v50 - v11 * v14 * v21 * v50
            - v12 * v16 * v21 * v50 - v12 * v24 * v50 + 2 * v1 * v2 * v21 * v4 * v5 * v50
            + 2 * v16 * v2 * v20 * v51 - 2 * v1 * v2 * v23 * v3 * v51
            + 4 * v0 * v11 * v17 * v4 * v51 + 2 * v13 * v23 * v4 * v51
            - 2 * v1 * v20 * v4 * v5 * v51 + 4 * v1 * v10 * v17 * v2 * v52
            + 2 * v15 * v2 * v21 * v52 - 4 * v0 * v1 * v20 * v4 * v52
            - 2 * v11 * v21 * v3 * v4 * v52 - 2 * v0 * v2 * v20 * v5 * v52
            - 2 * v1 * v2 * v21 * v3 * v5 * v52 + 2 * v13 * v21 * v4 * v5 * v52
            + 2 * v22 * v4 * v5 * v52 + 2 * v16 * v2 * v21 * v3 * v53 + 2 * v2 * v24 * v3 * v53
            + 4 * v1 * v10 * v17 * v4 * v53 + 2 * v15 * v21 * v4 * v53
            - 2 * v11 * v2 * v21 * v5 * v53 - 2 * v0 * v20 * v4 * v5 * v53
            - 2 * v1 * v21 * v3 * v4 * v5 * v53 - 4 * v10 * v11 * v17 * v54 - v13 * v16 * v21 * v54
            - v19 * v21 * v54 - v16 * v22 * v54 - v13 * v24 * v54 + 4 * v0 * v1 * v20 * v5 * v54
            + 2 * v11 * v21 * v3 * v5 * v54 + 2 * v0 * v12 * v46 * v55
            + 2 * v0 * v2 * v4 * v48 * v55 + 2 * v0 * v3 * v4 * v51 * v55 + 2 * v17 * v4 * v52 * v55
            + 2 * v10 * v2 * v53 * v55 + 2 * v1 * v2 * v4 * v45 * v56 + 2 * v1 * v18 * v49 * v56
            + 2 * v17 * v4 * v52 * v56 - 2 * v54 * v55 * v56 + 2 * v1 * v2 * v4 * v45 * v57
            + 2 * v0 * v12 * v46 * v57 + 2 * v0 * v2 * v4 * v48 * v57 + 2 * v1 * v18 * v49 * v57
            - 6 * v0 * v1 * v2 * v51 * v57 + 2 * v0 * v3 * v4 * v51 * v57 + 2 * v10 * v2 * v53 * v57
            - 2 * v54 * v55 * v57 - 2 * v54 * v56 * v57 - v11 * v45 * v58 - v17 * v47 * v58
            - v11 * v45 * v59 - 2 * v0 * v3 * v46 * v59 - 2 * v0 * v1 * v48 * v59 - v10 * v50 * v59
            + 2 * v0 * v1 * v14 * v3 * v45 * v6 + 2 * v15 * v18 * v46 * v6
            + 2 * v0 * v1 * v18 * v3 * v47 * v6 + 2 * v16 * v18 * v3 * v48 * v6
            + 2 * v12 * v23 * v49 * v6 + 2 * v0 * v1 * v12 * v45 * v5 * v6
            - 2 * v1 * v18 * v3 * v46 * v5 * v6 - 2 * v0 * v18 * v3 * v49 * v5 * v6
            + 2 * v2 * v23 * v4 * v50 * v6 + 2 * v0 * v1 * v18 * v5 * v50 * v6
            + 4 * v10 * v11 * v2 * v51 * v6 - 2 * v2 * v23 * v3 * v52 * v6
            - 4 * v1 * v2 * v23 * v53 * v6 - 2 * v23 * v3 * v4 * v53 * v6
            + 4 * v1 * v23 * v3 * v54 * v6 + 2 * v2 * v4 * v46 * v55 * v6 + 2 * v14 * v48 * v55 * v6
            + 2 * v2 * v5 * v51 * v55 * v6 - 6 * v0 * v2 * v52 * v55 * v6
            - 6 * v0 * v4 * v53 * v55 * v6 + 2 * v2 * v4 * v46 * v56 * v6 + 2 * v14 * v48 * v56 * v6
            - 6 * v1 * v4 * v51 * v56 * v6 + 2 * v2 * v5 * v51 * v56 * v6 - 2 * v1 * v46 * v58 * v6
            - 2 * v48 * v5 * v58 * v6 - 2 * v0 * v3 * v46 * v60 - v17 * v47 * v60
            + 2 * v3 * v45 * v5 * v60 - v10 * v50 * v60 - 2 * v0 * v49 * v6 * v60
            - 2 * v48 * v5 * v6 * v60 + v11 * v45 * v61 + v3 * v45 * v5 * v61 + v17 * v47 * v62
            + v3 * v45 * v5 * v62 + v3 * v45 * v5 * v63 + v10 * v50 * v63
            + 2 * v0 * v15 * v2 * v45 * v7 - 2 * v0 * v11 * v3 * v4 * v45 * v7
            - 2 * v10 * v13 * v4 * v46 * v7 - 4 * v11 * v17 * v4 * v46 * v7
            + 2 * v1 * v20 * v4 * v47 * v7 - 2 * v0 * v17 * v3 * v4 * v47 * v7
            - 4 * v10 * v11 * v2 * v48 * v7 - 2 * v16 * v17 * v2 * v48 * v7
            - 2 * v2 * v23 * v3 * v49 * v7 + 2 * v0 * v13 * v4 * v45 * v5 * v7
            - 2 * v20 * v4 * v49 * v5 * v7 + 2 * v1 * v2 * v23 * v50 * v7
            - 2 * v16 * v17 * v3 * v51 * v7 - 2 * v10 * v13 * v5 * v51 * v7
            + 4 * v0 * v11 * v17 * v52 * v7 + 2 * v13 * v23 * v52 * v7
            - 2 * v1 * v20 * v5 * v52 * v7 + 2 * v16 * v20 * v53 * v7 - 2 * v1 * v23 * v3 * v53 * v7
            + 2 * v3 * v4 * v46 * v55 * v7 + 2 * v0 * v2 * v49 * v55 * v7
            + 2 * v2 * v48 * v5 * v55 * v7 + 2 * v0 * v4 * v50 * v55 * v7
            - 4 * v3 * v5 * v51 * v55 * v7 + 2 * v0 * v3 * v52 * v55 * v7
            + 2 * v1 * v2 * v46 * v56 * v7 + 2 * v1 * v4 * v48 * v56 * v7
            + 2 * v2 * v48 * v5 * v56 * v7 + 2 * v11 * v51 * v56 * v7 + 2 * v1 * v2 * v46 * v57 * v7
            + 2 * v3 * v4 * v46 * v57 * v7 + 2 * v1 * v4 * v48 * v57 * v7
            + 2 * v0 * v2 * v49 * v57 * v7 + 2 * v0 * v4 * v50 * v57 * v7 + 2 * v11 * v51 * v57 * v7
            + 2 * v0 * v3 * v52 * v57 * v7 - 6 * v0 * v1 * v53 * v57 * v7
            + 2 * v16 * v2 * v3 * v45 * v6 * v7 + 2 * v15 * v4 * v45 * v6 * v7
            + 4 * v0 * v1 * v3 * v4 * v46 * v6 * v7 - 2 * v11 * v2 * v45 * v5 * v6 * v7
            - 4 * v0 * v2 * v3 * v46 * v5 * v6 * v7 + 4 * v0 * v1 * v2 * v48 * v5 * v6 * v7
            - 4 * v0 * v3 * v4 * v48 * v5 * v6 * v7 - 2 * v10 * v2 * v5 * v50 * v6 * v7
            - 4 * v0 * v15 * v51 * v6 * v7 + 4 * v0 * v1 * v3 * v5 * v51 * v6 * v7
            + 4 * v10 * v11 * v53 * v6 * v7 + 2 * v2 * v47 * v55 * v6 * v7
            + 2 * v4 * v49 * v55 * v6 * v7 + 2 * v5 * v53 * v55 * v6 * v7
            + 2 * v2 * v47 * v56 * v6 * v7 + 2 * v4 * v49 * v56 * v6 * v7
            - 6 * v1 * v52 * v56 * v6 * v7 + 2 * v5 * v53 * v56 * v6 * v7
            + 2 * v0 * v13 * v14 * v46 * v8 - 4 * v0 * v1 * v2 * v3 * v4 * v46 * v8
            - 2 * v1 * v17 * v2 * v4 * v47 * v8 - 2 * v0 * v1 * v14 * v3 * v48 * v8
            - 2 * v0 * v1 * v12 * v48 * v5 * v8 + 2 * v0 * v2 * v3 * v4 * v48 * v5 * v8
            - 2 * v1 * v10 * v2 * v4 * v50 * v8 - 2 * v0 * v15 * v2 * v51 * v8
            + 2 * v0 * v1 * v2 * v3 * v5 * v51 * v8 - 2 * v0 * v13 * v4 * v5 * v51 * v8
            - 2 * v10 * v13 * v4 * v52 * v8 - 2 * v16 * v17 * v2 * v53 * v8
            + 2 * v16 * v17 * v3 * v54 * v8 + 2 * v10 * v13 * v5 * v54 * v8
            + 2 * v2 * v4 * v49 * v56 * v8 + 2 * v1 * v2 * v52 * v56 * v8
            + 2 * v1 * v4 * v53 * v56 * v8 - 4 * v2 * v5 * v53 * v56 * v8
            + 2 * v2 * v4 * v49 * v57 * v8 + 2 * v1 * v2 * v52 * v57 * v8
            - 4 * v3 * v4 * v52 * v57 * v8 + 2 * v1 * v4 * v53 * v57 * v8 - 2 * v1 * v49 * v58 * v8
            + 2 * v5 * v50 * v58 * v8 + 2 * v3 * v47 * v59 * v8 - 2 * v1 * v49 * v59 * v8
            - 2 * v1 * v14 * v3 * v46 * v6 * v8 + 2 * v0 * v1 * v12 * v47 * v6 * v8
            + 2 * v12 * v16 * v48 * v6 * v8 - 2 * v0 * v14 * v3 * v49 * v6 * v8
            + 4 * v0 * v1 * v2 * v4 * v49 * v6 * v8 - 2 * v1 * v12 * v46 * v5 * v6 * v8
            + 2 * v2 * v3 * v4 * v46 * v5 * v6 * v8 - 4 * v1 * v2 * v4 * v48 * v5 * v6 * v8
            - 2 * v0 * v12 * v49 * v5 * v6 * v8 + 2 * v0 * v1 * v14 * v50 * v6 * v8
            - 2 * v16 * v2 * v3 * v51 * v6 * v8 - 2 * v15 * v4 * v51 * v6 * v8
            + 2 * v1 * v3 * v4 * v5 * v51 * v6 * v8 + 4 * v0 * v1 * v3 * v4 * v52 * v6 * v8
            + 2 * v0 * v2 * v3 * v5 * v52 * v6 * v8 + 4 * v0 * v1 * v2 * v5 * v53 * v6 * v8
            + 2 * v0 * v3 * v4 * v5 * v53 * v6 * v8 + 4 * v0 * v15 * v54 * v6 * v8
            - 4 * v0 * v1 * v3 * v5 * v54 * v6 * v8 + v3 * v47 * v61 * v8 + v5 * v50 * v61 * v8
            + v3 * v47 * v62 * v8 + v5 * v50 * v62 * v8 + v3 * v47 * v63 * v8 + v5 * v50 * v63 * v8
            - 2 * v15 * v2 * v46 * v7 * v8 + 2 * v0 * v13 * v4 * v47 * v7 * v8
            - 2 * v16 * v2 * v3 * v48 * v7 * v8 - 2 * v15 * v4 * v48 * v7 * v8
            - 4 * v0 * v1 * v3 * v4 * v49 * v7 * v8 + 2 * v1 * v2 * v3 * v46 * v5 * v7 * v8
            - 2 * v13 * v4 * v46 * v5 * v7 * v8 + 2 * v1 * v3 * v4 * v48 * v5 * v7 * v8
            + 2 * v0 * v2 * v3 * v49 * v5 * v7 * v8 + 2 * v13 * v16 * v51 * v7 * v8
            + 2 * v19 * v51 * v7 * v8 - 2 * v0 * v13 * v5 * v52 * v7 * v8
            - 2 * v0 * v15 * v53 * v7 * v8 + 2 * v0 * v1 * v3 * v5 * v53 * v7 * v8
            - 4 * v1 * v2 * v49 * v5 * v6 * v7 * v8 + 2 * v3 * v4 * v49 * v5 * v6 * v7 * v8
            + 2 * v16 * v2 * v50 * v6 * v7 * v8 - 2 * v15 * v52 * v6 * v7 * v8
            + 2 * v1 * v3 * v5 * v52 * v6 * v7 * v8 - 2 * v16 * v3 * v53 * v6 * v7 * v8))
        / 2.;
  gradalpha[1]
      = (v65 * v66
         * (-(v10 * v13 * v14 * v35) - v12 * v16 * v17 * v35 - v13 * v16 * v18 * v35
            - v18 * v19 * v35 + 2 * v1 * v14 * v20 * v36 - v10 * v13 * v18 * v37
            - v11 * v12 * v21 * v37 - v13 * v14 * v21 * v37 - v14 * v22 * v37
            + 2 * v0 * v15 * v18 * v38 + 2 * v1 * v12 * v23 * v38 + 2 * v0 * v14 * v20 * v39
            + 2 * v1 * v14 * v21 * v3 * v39 + 2 * v0 * v2 * v20 * v37 * v4
            + 2 * v1 * v2 * v21 * v3 * v37 * v4 - 2 * v2 * v23 * v3 * v38 * v4
            - 4 * v10 * v17 * v2 * v39 * v4 - 2 * v11 * v2 * v21 * v39 * v4 - v16 * v17 * v18 * v40
            - v11 * v14 * v21 * v40 - v12 * v16 * v21 * v40 - v12 * v24 * v40
            + 2 * v16 * v2 * v20 * v41 - 2 * v1 * v2 * v23 * v3 * v41
            + 4 * v0 * v11 * v17 * v4 * v41 + 2 * v13 * v23 * v4 * v41
            + 4 * v1 * v10 * v17 * v2 * v42 + 2 * v15 * v2 * v21 * v42
            - 4 * v0 * v1 * v20 * v4 * v42 - 2 * v11 * v21 * v3 * v4 * v42
            + 2 * v16 * v2 * v21 * v3 * v43 + 2 * v2 * v24 * v3 * v43
            + 4 * v1 * v10 * v17 * v4 * v43 + 2 * v15 * v21 * v4 * v43 - 4 * v10 * v11 * v17 * v44
            - v13 * v16 * v21 * v44 - v19 * v21 * v44 - v16 * v22 * v44 - v13 * v24 * v44
            + 2 * v0 * v13 * v18 * v36 * v5 - 2 * v0 * v1 * v18 * v3 * v38 * v5
            + 2 * v1 * v12 * v21 * v39 * v5 - 2 * v2 * v20 * v36 * v4 * v5
            - 2 * v2 * v21 * v3 * v39 * v4 * v5 + 2 * v1 * v2 * v21 * v4 * v40 * v5
            - 2 * v1 * v20 * v4 * v41 * v5 - 2 * v0 * v2 * v20 * v42 * v5
            - 2 * v1 * v2 * v21 * v3 * v42 * v5 + 2 * v13 * v21 * v4 * v42 * v5
            + 2 * v22 * v4 * v42 * v5 - 2 * v11 * v2 * v21 * v43 * v5 - 2 * v0 * v20 * v4 * v43 * v5
            - 2 * v1 * v21 * v3 * v4 * v43 * v5 + 4 * v0 * v1 * v20 * v44 * v5
            + 2 * v11 * v21 * v3 * v44 * v5 + 2 * v0 * v12 * v36 * v55
            + 2 * v0 * v2 * v38 * v4 * v55 + 2 * v0 * v3 * v4 * v41 * v55 + 2 * v17 * v4 * v42 * v55
            + 2 * v10 * v2 * v43 * v55 + 2 * v1 * v18 * v39 * v56 + 2 * v1 * v2 * v35 * v4 * v56
            + 2 * v17 * v4 * v42 * v56 - 2 * v44 * v55 * v56 + 2 * v0 * v12 * v36 * v57
            + 2 * v1 * v18 * v39 * v57 + 2 * v1 * v2 * v35 * v4 * v57 + 2 * v0 * v2 * v38 * v4 * v57
            - 6 * v0 * v1 * v2 * v41 * v57 + 2 * v0 * v3 * v4 * v41 * v57 + 2 * v10 * v2 * v43 * v57
            - 2 * v44 * v55 * v57 - 2 * v44 * v56 * v57 - v11 * v35 * v58 - v17 * v37 * v58
            - v11 * v35 * v59 - 2 * v0 * v3 * v36 * v59 - 2 * v0 * v1 * v38 * v59 - v10 * v40 * v59
            + 2 * v0 * v1 * v14 * v3 * v35 * v6 + 2 * v15 * v18 * v36 * v6
            + 2 * v0 * v1 * v18 * v3 * v37 * v6 + 2 * v16 * v18 * v3 * v38 * v6
            + 2 * v12 * v23 * v39 * v6 + 2 * v2 * v23 * v4 * v40 * v6
            + 4 * v10 * v11 * v2 * v41 * v6 - 2 * v2 * v23 * v3 * v42 * v6
            - 4 * v1 * v2 * v23 * v43 * v6 - 2 * v23 * v3 * v4 * v43 * v6
            + 4 * v1 * v23 * v3 * v44 * v6 + 2 * v0 * v1 * v12 * v35 * v5 * v6
            - 2 * v1 * v18 * v3 * v36 * v5 * v6 - 2 * v0 * v18 * v3 * v39 * v5 * v6
            + 2 * v0 * v1 * v18 * v40 * v5 * v6 + 2 * v14 * v38 * v55 * v6
            + 2 * v2 * v36 * v4 * v55 * v6 - 6 * v0 * v2 * v42 * v55 * v6
            - 6 * v0 * v4 * v43 * v55 * v6 + 2 * v2 * v41 * v5 * v55 * v6 + 2 * v14 * v38 * v56 * v6
            + 2 * v2 * v36 * v4 * v56 * v6 - 6 * v1 * v4 * v41 * v56 * v6
            + 2 * v2 * v41 * v5 * v56 * v6 - 2 * v1 * v36 * v58 * v6 - 2 * v38 * v5 * v58 * v6
            - 2 * v0 * v3 * v36 * v60 - v17 * v37 * v60 - v10 * v40 * v60 + 2 * v3 * v35 * v5 * v60
            - 2 * v0 * v39 * v6 * v60 - 2 * v38 * v5 * v6 * v60 + v11 * v35 * v61
            + v3 * v35 * v5 * v61 + v17 * v37 * v62 + v3 * v35 * v5 * v62 + v10 * v40 * v63
            + v3 * v35 * v5 * v63 + 2 * v0 * v15 * v2 * v35 * v7 - 4 * v10 * v11 * v2 * v38 * v7
            - 2 * v16 * v17 * v2 * v38 * v7 - 2 * v2 * v23 * v3 * v39 * v7
            - 2 * v0 * v11 * v3 * v35 * v4 * v7 - 2 * v10 * v13 * v36 * v4 * v7
            - 4 * v11 * v17 * v36 * v4 * v7 + 2 * v1 * v20 * v37 * v4 * v7
            - 2 * v0 * v17 * v3 * v37 * v4 * v7 + 2 * v1 * v2 * v23 * v40 * v7
            - 2 * v16 * v17 * v3 * v41 * v7 + 4 * v0 * v11 * v17 * v42 * v7
            + 2 * v13 * v23 * v42 * v7 + 2 * v16 * v20 * v43 * v7 - 2 * v1 * v23 * v3 * v43 * v7
            + 2 * v0 * v13 * v35 * v4 * v5 * v7 - 2 * v20 * v39 * v4 * v5 * v7
            - 2 * v10 * v13 * v41 * v5 * v7 - 2 * v1 * v20 * v42 * v5 * v7
            + 2 * v0 * v2 * v39 * v55 * v7 + 2 * v3 * v36 * v4 * v55 * v7
            + 2 * v0 * v4 * v40 * v55 * v7 + 2 * v0 * v3 * v42 * v55 * v7
            + 2 * v2 * v38 * v5 * v55 * v7 - 4 * v3 * v41 * v5 * v55 * v7
            + 2 * v1 * v2 * v36 * v56 * v7 + 2 * v1 * v38 * v4 * v56 * v7 + 2 * v11 * v41 * v56 * v7
            + 2 * v2 * v38 * v5 * v56 * v7 + 2 * v1 * v2 * v36 * v57 * v7
            + 2 * v0 * v2 * v39 * v57 * v7 + 2 * v3 * v36 * v4 * v57 * v7
            + 2 * v1 * v38 * v4 * v57 * v7 + 2 * v0 * v4 * v40 * v57 * v7 + 2 * v11 * v41 * v57 * v7
            + 2 * v0 * v3 * v42 * v57 * v7 - 6 * v0 * v1 * v43 * v57 * v7
            + 2 * v16 * v2 * v3 * v35 * v6 * v7 + 2 * v15 * v35 * v4 * v6 * v7
            + 4 * v0 * v1 * v3 * v36 * v4 * v6 * v7 - 4 * v0 * v15 * v41 * v6 * v7
            + 4 * v10 * v11 * v43 * v6 * v7 - 2 * v11 * v2 * v35 * v5 * v6 * v7
            - 4 * v0 * v2 * v3 * v36 * v5 * v6 * v7 + 4 * v0 * v1 * v2 * v38 * v5 * v6 * v7
            - 4 * v0 * v3 * v38 * v4 * v5 * v6 * v7 - 2 * v10 * v2 * v40 * v5 * v6 * v7
            + 4 * v0 * v1 * v3 * v41 * v5 * v6 * v7 + 2 * v2 * v37 * v55 * v6 * v7
            + 2 * v39 * v4 * v55 * v6 * v7 + 2 * v43 * v5 * v55 * v6 * v7
            + 2 * v2 * v37 * v56 * v6 * v7 + 2 * v39 * v4 * v56 * v6 * v7
            - 6 * v1 * v42 * v56 * v6 * v7 + 2 * v43 * v5 * v56 * v6 * v7
            + 2 * v0 * v13 * v14 * v36 * v8 - 2 * v0 * v1 * v14 * v3 * v38 * v8
            - 4 * v0 * v1 * v2 * v3 * v36 * v4 * v8 - 2 * v1 * v17 * v2 * v37 * v4 * v8
            - 2 * v1 * v10 * v2 * v4 * v40 * v8 - 2 * v0 * v15 * v2 * v41 * v8
            - 2 * v10 * v13 * v4 * v42 * v8 - 2 * v16 * v17 * v2 * v43 * v8
            + 2 * v16 * v17 * v3 * v44 * v8 - 2 * v0 * v1 * v12 * v38 * v5 * v8
            + 2 * v0 * v2 * v3 * v38 * v4 * v5 * v8 + 2 * v0 * v1 * v2 * v3 * v41 * v5 * v8
            - 2 * v0 * v13 * v4 * v41 * v5 * v8 + 2 * v10 * v13 * v44 * v5 * v8
            + 2 * v2 * v39 * v4 * v56 * v8 + 2 * v1 * v2 * v42 * v56 * v8
            + 2 * v1 * v4 * v43 * v56 * v8 - 4 * v2 * v43 * v5 * v56 * v8
            + 2 * v2 * v39 * v4 * v57 * v8 + 2 * v1 * v2 * v42 * v57 * v8
            - 4 * v3 * v4 * v42 * v57 * v8 + 2 * v1 * v4 * v43 * v57 * v8 - 2 * v1 * v39 * v58 * v8
            + 2 * v40 * v5 * v58 * v8 + 2 * v3 * v37 * v59 * v8 - 2 * v1 * v39 * v59 * v8
            - 2 * v1 * v14 * v3 * v36 * v6 * v8 + 2 * v0 * v1 * v12 * v37 * v6 * v8
            + 2 * v12 * v16 * v38 * v6 * v8 - 2 * v0 * v14 * v3 * v39 * v6 * v8
            + 4 * v0 * v1 * v2 * v39 * v4 * v6 * v8 + 2 * v0 * v1 * v14 * v40 * v6 * v8
            - 2 * v16 * v2 * v3 * v41 * v6 * v8 - 2 * v15 * v4 * v41 * v6 * v8
            + 4 * v0 * v1 * v3 * v4 * v42 * v6 * v8 + 4 * v0 * v15 * v44 * v6 * v8
            - 2 * v1 * v12 * v36 * v5 * v6 * v8 - 2 * v0 * v12 * v39 * v5 * v6 * v8
            + 2 * v2 * v3 * v36 * v4 * v5 * v6 * v8 - 4 * v1 * v2 * v38 * v4 * v5 * v6 * v8
            + 2 * v1 * v3 * v4 * v41 * v5 * v6 * v8 + 2 * v0 * v2 * v3 * v42 * v5 * v6 * v8
            + 4 * v0 * v1 * v2 * v43 * v5 * v6 * v8 + 2 * v0 * v3 * v4 * v43 * v5 * v6 * v8
            - 4 * v0 * v1 * v3 * v44 * v5 * v6 * v8 + v3 * v37 * v61 * v8 + v40 * v5 * v61 * v8
            + v3 * v37 * v62 * v8 + v40 * v5 * v62 * v8 + v3 * v37 * v63 * v8 + v40 * v5 * v63 * v8
            - 2 * v15 * v2 * v36 * v7 * v8 - 2 * v16 * v2 * v3 * v38 * v7 * v8
            + 2 * v0 * v13 * v37 * v4 * v7 * v8 - 2 * v15 * v38 * v4 * v7 * v8
            - 4 * v0 * v1 * v3 * v39 * v4 * v7 * v8 + 2 * v13 * v16 * v41 * v7 * v8
            + 2 * v19 * v41 * v7 * v8 - 2 * v0 * v15 * v43 * v7 * v8
            + 2 * v1 * v2 * v3 * v36 * v5 * v7 * v8 + 2 * v0 * v2 * v3 * v39 * v5 * v7 * v8
            - 2 * v13 * v36 * v4 * v5 * v7 * v8 + 2 * v1 * v3 * v38 * v4 * v5 * v7 * v8
            - 2 * v0 * v13 * v42 * v5 * v7 * v8 + 2 * v0 * v1 * v3 * v43 * v5 * v7 * v8
            + 2 * v16 * v2 * v40 * v6 * v7 * v8 - 2 * v15 * v42 * v6 * v7 * v8
            - 2 * v16 * v3 * v43 * v6 * v7 * v8 - 4 * v1 * v2 * v39 * v5 * v6 * v7 * v8
            + 2 * v3 * v39 * v4 * v5 * v6 * v7 * v8 + 2 * v1 * v3 * v42 * v5 * v6 * v7 * v8))
        / 2.;
  gradalpha[2]
      = (v65 * v66
         * (-(v10 * v13 * v14 * v25) - v12 * v16 * v17 * v25 - v13 * v16 * v18 * v25
            - v18 * v19 * v25 + 2 * v1 * v14 * v20 * v26 - v10 * v13 * v18 * v27
            - v11 * v12 * v21 * v27 - v13 * v14 * v21 * v27 - v14 * v22 * v27
            + 2 * v0 * v15 * v18 * v28 + 2 * v1 * v12 * v23 * v28 + 2 * v0 * v14 * v20 * v29
            + 2 * v1 * v14 * v21 * v29 * v3 - v16 * v17 * v18 * v30 - v11 * v14 * v21 * v30
            - v12 * v16 * v21 * v30 - v12 * v24 * v30 + 2 * v16 * v2 * v20 * v31
            - 2 * v1 * v2 * v23 * v3 * v31 + 4 * v1 * v10 * v17 * v2 * v32
            + 2 * v15 * v2 * v21 * v32 + 2 * v16 * v2 * v21 * v3 * v33 + 2 * v2 * v24 * v3 * v33
            - 4 * v10 * v11 * v17 * v34 - v13 * v16 * v21 * v34 - v19 * v21 * v34 - v16 * v22 * v34
            - v13 * v24 * v34 + 2 * v0 * v2 * v20 * v27 * v4 - 4 * v10 * v17 * v2 * v29 * v4
            - 2 * v11 * v2 * v21 * v29 * v4 + 2 * v1 * v2 * v21 * v27 * v3 * v4
            - 2 * v2 * v23 * v28 * v3 * v4 + 4 * v0 * v11 * v17 * v31 * v4
            + 2 * v13 * v23 * v31 * v4 - 4 * v0 * v1 * v20 * v32 * v4
            - 2 * v11 * v21 * v3 * v32 * v4 + 4 * v1 * v10 * v17 * v33 * v4
            + 2 * v15 * v21 * v33 * v4 + 2 * v0 * v13 * v18 * v26 * v5
            + 2 * v1 * v12 * v21 * v29 * v5 - 2 * v0 * v1 * v18 * v28 * v3 * v5
            - 2 * v0 * v2 * v20 * v32 * v5 - 2 * v1 * v2 * v21 * v3 * v32 * v5
            - 2 * v11 * v2 * v21 * v33 * v5 + 4 * v0 * v1 * v20 * v34 * v5
            + 2 * v11 * v21 * v3 * v34 * v5 - 2 * v2 * v20 * v26 * v4 * v5
            - 2 * v2 * v21 * v29 * v3 * v4 * v5 + 2 * v1 * v2 * v21 * v30 * v4 * v5
            - 2 * v1 * v20 * v31 * v4 * v5 + 2 * v13 * v21 * v32 * v4 * v5 + 2 * v22 * v32 * v4 * v5
            - 2 * v0 * v20 * v33 * v4 * v5 - 2 * v1 * v21 * v3 * v33 * v4 * v5
            + 2 * v0 * v12 * v26 * v55 + 2 * v10 * v2 * v33 * v55 + 2 * v0 * v2 * v28 * v4 * v55
            + 2 * v0 * v3 * v31 * v4 * v55 + 2 * v17 * v32 * v4 * v55 + 2 * v1 * v18 * v29 * v56
            + 2 * v1 * v2 * v25 * v4 * v56 + 2 * v17 * v32 * v4 * v56 - 2 * v34 * v55 * v56
            + 2 * v0 * v12 * v26 * v57 + 2 * v1 * v18 * v29 * v57 - 6 * v0 * v1 * v2 * v31 * v57
            + 2 * v10 * v2 * v33 * v57 + 2 * v1 * v2 * v25 * v4 * v57 + 2 * v0 * v2 * v28 * v4 * v57
            + 2 * v0 * v3 * v31 * v4 * v57 - 2 * v34 * v55 * v57 - 2 * v34 * v56 * v57
            - v11 * v25 * v58 - v17 * v27 * v58 - v11 * v25 * v59 - 2 * v0 * v1 * v28 * v59
            - 2 * v0 * v26 * v3 * v59 - v10 * v30 * v59 + 2 * v15 * v18 * v26 * v6
            + 2 * v12 * v23 * v29 * v6 + 2 * v0 * v1 * v14 * v25 * v3 * v6
            + 2 * v0 * v1 * v18 * v27 * v3 * v6 + 2 * v16 * v18 * v28 * v3 * v6
            + 4 * v10 * v11 * v2 * v31 * v6 - 2 * v2 * v23 * v3 * v32 * v6
            - 4 * v1 * v2 * v23 * v33 * v6 + 4 * v1 * v23 * v3 * v34 * v6
            + 2 * v2 * v23 * v30 * v4 * v6 - 2 * v23 * v3 * v33 * v4 * v6
            + 2 * v0 * v1 * v12 * v25 * v5 * v6 - 2 * v1 * v18 * v26 * v3 * v5 * v6
            - 2 * v0 * v18 * v29 * v3 * v5 * v6 + 2 * v0 * v1 * v18 * v30 * v5 * v6
            + 2 * v14 * v28 * v55 * v6 - 6 * v0 * v2 * v32 * v55 * v6 + 2 * v2 * v26 * v4 * v55 * v6
            - 6 * v0 * v33 * v4 * v55 * v6 + 2 * v2 * v31 * v5 * v55 * v6 + 2 * v14 * v28 * v56 * v6
            + 2 * v2 * v26 * v4 * v56 * v6 - 6 * v1 * v31 * v4 * v56 * v6
            + 2 * v2 * v31 * v5 * v56 * v6 - 2 * v1 * v26 * v58 * v6 - 2 * v28 * v5 * v58 * v6
            - v17 * v27 * v60 - 2 * v0 * v26 * v3 * v60 - v10 * v30 * v60 + 2 * v25 * v3 * v5 * v60
            - 2 * v0 * v29 * v6 * v60 - 2 * v28 * v5 * v6 * v60 + v11 * v25 * v61
            + v25 * v3 * v5 * v61 + v17 * v27 * v62 + v25 * v3 * v5 * v62 + v10 * v30 * v63
            + v25 * v3 * v5 * v63 + 2 * v0 * v15 * v2 * v25 * v7 - 4 * v10 * v11 * v2 * v28 * v7
            - 2 * v16 * v17 * v2 * v28 * v7 - 2 * v2 * v23 * v29 * v3 * v7
            + 2 * v1 * v2 * v23 * v30 * v7 - 2 * v16 * v17 * v3 * v31 * v7
            + 4 * v0 * v11 * v17 * v32 * v7 + 2 * v13 * v23 * v32 * v7 + 2 * v16 * v20 * v33 * v7
            - 2 * v1 * v23 * v3 * v33 * v7 - 2 * v10 * v13 * v26 * v4 * v7
            - 4 * v11 * v17 * v26 * v4 * v7 + 2 * v1 * v20 * v27 * v4 * v7
            - 2 * v0 * v11 * v25 * v3 * v4 * v7 - 2 * v0 * v17 * v27 * v3 * v4 * v7
            - 2 * v10 * v13 * v31 * v5 * v7 - 2 * v1 * v20 * v32 * v5 * v7
            + 2 * v0 * v13 * v25 * v4 * v5 * v7 - 2 * v20 * v29 * v4 * v5 * v7
            + 2 * v0 * v2 * v29 * v55 * v7 + 2 * v0 * v3 * v32 * v55 * v7
            + 2 * v26 * v3 * v4 * v55 * v7 + 2 * v0 * v30 * v4 * v55 * v7
            + 2 * v2 * v28 * v5 * v55 * v7 - 4 * v3 * v31 * v5 * v55 * v7
            + 2 * v1 * v2 * v26 * v56 * v7 + 2 * v11 * v31 * v56 * v7 + 2 * v1 * v28 * v4 * v56 * v7
            + 2 * v2 * v28 * v5 * v56 * v7 + 2 * v1 * v2 * v26 * v57 * v7
            + 2 * v0 * v2 * v29 * v57 * v7 + 2 * v11 * v31 * v57 * v7 + 2 * v0 * v3 * v32 * v57 * v7
            - 6 * v0 * v1 * v33 * v57 * v7 + 2 * v1 * v28 * v4 * v57 * v7
            + 2 * v26 * v3 * v4 * v57 * v7 + 2 * v0 * v30 * v4 * v57 * v7
            + 2 * v16 * v2 * v25 * v3 * v6 * v7 - 4 * v0 * v15 * v31 * v6 * v7
            + 4 * v10 * v11 * v33 * v6 * v7 + 2 * v15 * v25 * v4 * v6 * v7
            + 4 * v0 * v1 * v26 * v3 * v4 * v6 * v7 - 2 * v11 * v2 * v25 * v5 * v6 * v7
            + 4 * v0 * v1 * v2 * v28 * v5 * v6 * v7 - 4 * v0 * v2 * v26 * v3 * v5 * v6 * v7
            - 2 * v10 * v2 * v30 * v5 * v6 * v7 + 4 * v0 * v1 * v3 * v31 * v5 * v6 * v7
            - 4 * v0 * v28 * v3 * v4 * v5 * v6 * v7 + 2 * v2 * v27 * v55 * v6 * v7
            + 2 * v29 * v4 * v55 * v6 * v7 + 2 * v33 * v5 * v55 * v6 * v7
            + 2 * v2 * v27 * v56 * v6 * v7 - 6 * v1 * v32 * v56 * v6 * v7
            + 2 * v29 * v4 * v56 * v6 * v7 + 2 * v33 * v5 * v56 * v6 * v7
            + 2 * v0 * v13 * v14 * v26 * v8 - 2 * v0 * v1 * v14 * v28 * v3 * v8
            - 2 * v0 * v15 * v2 * v31 * v8 - 2 * v16 * v17 * v2 * v33 * v8
            + 2 * v16 * v17 * v3 * v34 * v8 - 2 * v1 * v17 * v2 * v27 * v4 * v8
            - 4 * v0 * v1 * v2 * v26 * v3 * v4 * v8 - 2 * v1 * v10 * v2 * v30 * v4 * v8
            - 2 * v10 * v13 * v32 * v4 * v8 - 2 * v0 * v1 * v12 * v28 * v5 * v8
            + 2 * v0 * v1 * v2 * v3 * v31 * v5 * v8 + 2 * v10 * v13 * v34 * v5 * v8
            + 2 * v0 * v2 * v28 * v3 * v4 * v5 * v8 - 2 * v0 * v13 * v31 * v4 * v5 * v8
            + 2 * v1 * v2 * v32 * v56 * v8 + 2 * v2 * v29 * v4 * v56 * v8
            + 2 * v1 * v33 * v4 * v56 * v8 - 4 * v2 * v33 * v5 * v56 * v8
            + 2 * v1 * v2 * v32 * v57 * v8 + 2 * v2 * v29 * v4 * v57 * v8
            - 4 * v3 * v32 * v4 * v57 * v8 + 2 * v1 * v33 * v4 * v57 * v8 - 2 * v1 * v29 * v58 * v8
            + 2 * v30 * v5 * v58 * v8 - 2 * v1 * v29 * v59 * v8 + 2 * v27 * v3 * v59 * v8
            + 2 * v0 * v1 * v12 * v27 * v6 * v8 + 2 * v12 * v16 * v28 * v6 * v8
            - 2 * v1 * v14 * v26 * v3 * v6 * v8 - 2 * v0 * v14 * v29 * v3 * v6 * v8
            + 2 * v0 * v1 * v14 * v30 * v6 * v8 - 2 * v16 * v2 * v3 * v31 * v6 * v8
            + 4 * v0 * v15 * v34 * v6 * v8 + 4 * v0 * v1 * v2 * v29 * v4 * v6 * v8
            - 2 * v15 * v31 * v4 * v6 * v8 + 4 * v0 * v1 * v3 * v32 * v4 * v6 * v8
            - 2 * v1 * v12 * v26 * v5 * v6 * v8 - 2 * v0 * v12 * v29 * v5 * v6 * v8
            + 2 * v0 * v2 * v3 * v32 * v5 * v6 * v8 + 4 * v0 * v1 * v2 * v33 * v5 * v6 * v8
            - 4 * v0 * v1 * v3 * v34 * v5 * v6 * v8 - 4 * v1 * v2 * v28 * v4 * v5 * v6 * v8
            + 2 * v2 * v26 * v3 * v4 * v5 * v6 * v8 + 2 * v1 * v3 * v31 * v4 * v5 * v6 * v8
            + 2 * v0 * v3 * v33 * v4 * v5 * v6 * v8 + v27 * v3 * v61 * v8 + v30 * v5 * v61 * v8
            + v27 * v3 * v62 * v8 + v30 * v5 * v62 * v8 + v27 * v3 * v63 * v8 + v30 * v5 * v63 * v8
            - 2 * v15 * v2 * v26 * v7 * v8 - 2 * v16 * v2 * v28 * v3 * v7 * v8
            + 2 * v13 * v16 * v31 * v7 * v8 + 2 * v19 * v31 * v7 * v8 - 2 * v0 * v15 * v33 * v7 * v8
            + 2 * v0 * v13 * v27 * v4 * v7 * v8 - 2 * v15 * v28 * v4 * v7 * v8
            - 4 * v0 * v1 * v29 * v3 * v4 * v7 * v8 + 2 * v1 * v2 * v26 * v3 * v5 * v7 * v8
            + 2 * v0 * v2 * v29 * v3 * v5 * v7 * v8 - 2 * v0 * v13 * v32 * v5 * v7 * v8
            + 2 * v0 * v1 * v3 * v33 * v5 * v7 * v8 - 2 * v13 * v26 * v4 * v5 * v7 * v8
            + 2 * v1 * v28 * v3 * v4 * v5 * v7 * v8 + 2 * v16 * v2 * v30 * v6 * v7 * v8
            - 2 * v15 * v32 * v6 * v7 * v8 - 2 * v16 * v3 * v33 * v6 * v7 * v8
            - 4 * v1 * v2 * v29 * v5 * v6 * v7 * v8 + 2 * v1 * v3 * v32 * v5 * v6 * v7 * v8
            + 2 * v29 * v3 * v4 * v5 * v6 * v7 * v8))
        / 2.;

  return gradalpha;
}