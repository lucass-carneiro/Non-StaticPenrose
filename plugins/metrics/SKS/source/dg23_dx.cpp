#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::dllgSKS_23_dx(double t, double x, double y, double z) const noexcept
    -> double {

  using namespace sks_aux;
  using std::pow;
  using std::sqrt;

  const double v0{Z(2, t, x, y, z)};
  const double v1{Y(2, t, x, y, z)};
  const double v2{X(2, t, x, y, z)};
  const double v3{Z(1, t, x, y, z)};
  const double v4{Y(1, t, x, y, z)};
  const double v5{X(1, t, x, y, z)};
  const double v6{dsy_dt(2, t)};
  const double v7{dsx_dt(2, t)};
  const double v8{dsy_dt(1, t)};
  const double v9{dsx_dt(1, t)};
  const double v10{pow(v6, 2)};
  const double v11{pow(v7, 2)};
  const double v12{pow(v8, 2)};
  const double v13{pow(v9, 2)};
  const double v14{pow(v6, 3)};
  const double v15{pow(v8, 3)};
  const double v16{dX_dx(2, t, x, y, z)};
  const double v17{-v10};
  const double v18{-v11};
  const double v19{-v12};
  const double v20{-v13};
  const double v21{dX_dx(1, t, x, y, z)};
  const double v22{dZ_dx(2, t, x, y, z)};
  const double v23{dZ_dx(1, t, x, y, z)};
  const double v24{dY_dx(2, t, x, y, z)};
  const double v25{dY_dx(1, t, x, y, z)};
  const double v26{l3_KS(a2, v2, v1, v0)};
  const double v27{l1_KS(a2, v2, v1, v0)};
  const double v28{l2_KS(a2, v2, v1, v0)};
  const double v29{l3_KS(a1, v5, v4, v3)};
  const double v30{l1_KS(a1, v5, v4, v3)};
  const double v31{l2_KS(a1, v5, v4, v3)};
  const double v32{H_KS(M2, a2, v2, v1, v0)};
  const double v33{H_KS(M1, a1, v5, v4, v3)};
  const double v34{dl3_KS_dZ(a2, v2, v1, v0)};
  const double v35{dl2_KS_dZ(a2, v2, v1, v0)};
  const double v36{dl1_KS_dZ(a2, v2, v1, v0)};
  const double v37{dl3_KS_dZ(a1, v5, v4, v3)};
  const double v38{dl2_KS_dZ(a1, v5, v4, v3)};
  const double v39{dl1_KS_dZ(a1, v5, v4, v3)};
  const double v40{dl3_KS_dY(a2, v2, v1, v0)};
  const double v41{dl2_KS_dY(a2, v2, v1, v0)};
  const double v42{dl1_KS_dY(a2, v2, v1, v0)};
  const double v43{dl3_KS_dY(a1, v5, v4, v3)};
  const double v44{dl2_KS_dY(a1, v5, v4, v3)};
  const double v45{dl1_KS_dY(a1, v5, v4, v3)};
  const double v46{dl3_KS_dX(a2, v2, v1, v0)};
  const double v47{dl2_KS_dX(a2, v2, v1, v0)};
  const double v48{dl1_KS_dX(a2, v2, v1, v0)};
  const double v49{dl3_KS_dX(a1, v5, v4, v3)};
  const double v50{dl2_KS_dX(a1, v5, v4, v3)};
  const double v51{dl1_KS_dX(a1, v5, v4, v3)};
  const double v52{1 + v17 + v18};
  const double v53{1 + v19 + v20};
  const double v54{dH_KS_dX(M2, a2, v2, v1, v0)};
  const double v55{dH_KS_dX(M1, a1, v5, v4, v3)};
  const double v56{dH_KS_dZ(M2, a2, v2, v1, v0)};
  const double v57{dH_KS_dZ(M1, a1, v5, v4, v3)};
  const double v58{dH_KS_dY(M2, a2, v2, v1, v0)};
  const double v59{dH_KS_dY(M1, a1, v5, v4, v3)};
  const double v60{sqrt(v52)};
  const double v61{sqrt(v53)};

  return (v14 * sqrt(1 - v13 + v19) * v20 * v22 * v32 * v34
          + v14 * v19 * sqrt(1 - v12 + v20) * v22 * v32 * v34
          + v15 * sqrt(1 - v11 + v17) * v18 * v23 * v33 * v37
          + v15 * v17 * sqrt(1 - v10 + v18) * v23 * v33 * v37
          + v14 * sqrt(1 - v13 + v19) * v20 * v24 * v32 * v40
          + v14 * v19 * sqrt(1 - v12 + v20) * v24 * v32 * v40
          + v15 * sqrt(1 - v11 + v17) * v18 * v25 * v33 * v43
          + v15 * v17 * sqrt(1 - v10 + v18) * v25 * v33 * v43
          + v14 * v16 * sqrt(1 - v13 + v19) * v20 * v32 * v46
          + v14 * v16 * v19 * sqrt(1 - v12 + v20) * v32 * v46
          + v15 * sqrt(1 - v11 + v17) * v18 * v21 * v33 * v49
          + v15 * v17 * sqrt(1 - v10 + v18) * v21 * v33 * v49
          + v14 * v16 * sqrt(1 - v13 + v19) * v20 * v26 * v54
          + v14 * v16 * v19 * sqrt(1 - v12 + v20) * v26 * v54
          + v15 * sqrt(1 - v11 + v17) * v18 * v21 * v29 * v55
          + v15 * v17 * sqrt(1 - v10 + v18) * v21 * v29 * v55
          + v14 * sqrt(1 - v13 + v19) * v20 * v22 * v26 * v56
          + v14 * v19 * sqrt(1 - v12 + v20) * v22 * v26 * v56
          + v15 * sqrt(1 - v11 + v17) * v18 * v23 * v29 * v57
          + v15 * v17 * sqrt(1 - v10 + v18) * v23 * v29 * v57
          + v14 * sqrt(1 - v13 + v19) * v20 * v24 * v26 * v58
          + v14 * v19 * sqrt(1 - v12 + v20) * v24 * v26 * v58
          + v15 * sqrt(1 - v11 + v17) * v18 * v25 * v29 * v59
          + v15 * v17 * sqrt(1 - v10 + v18) * v25 * v29 * v59
          + v10 * v12 * v23 * v31 * v33 * v37 * v60 + v11 * v12 * v23 * v31 * v33 * v37 * v60
          + v10 * v12 * v23 * v29 * v33 * v38 * v60 + v11 * v12 * v23 * v29 * v33 * v38 * v60
          + v10 * v12 * v25 * v31 * v33 * v43 * v60 + v11 * v12 * v25 * v31 * v33 * v43 * v60
          + v10 * v12 * v25 * v29 * v33 * v44 * v60 + v11 * v12 * v25 * v29 * v33 * v44 * v60
          + v10 * v12 * v21 * v31 * v33 * v49 * v60 + v11 * v12 * v21 * v31 * v33 * v49 * v60
          + v10 * v12 * v21 * v29 * v33 * v50 * v60 + v11 * v12 * v21 * v29 * v33 * v50 * v60
          + v10 * v12 * v21 * v29 * v31 * v55 * v60 + v11 * v12 * v21 * v29 * v31 * v55 * v60
          + v10 * v12 * v23 * v29 * v31 * v57 * v60 + v11 * v12 * v23 * v29 * v31 * v57 * v60
          + v10 * v12 * v25 * v29 * v31 * v59 * v60 + v11 * v12 * v25 * v29 * v31 * v59 * v60
          + v10 * v12 * v22 * v28 * v32 * v34 * v61 + v10 * v13 * v22 * v28 * v32 * v34 * v61
          + v10 * v12 * v22 * v26 * v32 * v35 * v61 + v10 * v13 * v22 * v26 * v32 * v35 * v61
          + v10 * v12 * v24 * v28 * v32 * v40 * v61 + v10 * v13 * v24 * v28 * v32 * v40 * v61
          + v10 * v12 * v24 * v26 * v32 * v41 * v61 + v10 * v13 * v24 * v26 * v32 * v41 * v61
          + v10 * v12 * v16 * v28 * v32 * v46 * v61 + v10 * v13 * v16 * v28 * v32 * v46 * v61
          + v10 * v12 * v16 * v26 * v32 * v47 * v61 + v10 * v13 * v16 * v26 * v32 * v47 * v61
          + v10 * v12 * v16 * v26 * v28 * v54 * v61 + v10 * v13 * v16 * v26 * v28 * v54 * v61
          + v10 * v12 * v22 * v26 * v28 * v56 * v61 + v10 * v13 * v22 * v26 * v28 * v56 * v61
          + v10 * v12 * v24 * v26 * v28 * v58 * v61 + v10 * v13 * v24 * v26 * v28 * v58 * v61
          + v12 * v18 * v22 * v32 * v34 * v6 * v61 + v13 * v18 * v22 * v32 * v34 * v6 * v61
          + v12 * v18 * v24 * v32 * v40 * v6 * v61 + v13 * v18 * v24 * v32 * v40 * v6 * v61
          + v12 * v16 * v18 * v32 * v46 * v6 * v61 + v13 * v16 * v18 * v32 * v46 * v6 * v61
          + v12 * v16 * v18 * v26 * v54 * v6 * v61 + v13 * v16 * v18 * v26 * v54 * v6 * v61
          + v12 * v18 * v22 * v26 * v56 * v6 * v61 + v13 * v18 * v22 * v26 * v56 * v6 * v61
          + v12 * v18 * v24 * v26 * v58 * v6 * v61 + v13 * v18 * v24 * v26 * v58 * v6 * v61
          + v11 * v12 * v22 * v28 * v32 * v34 * v60 * v61
          + v11 * v13 * v22 * v28 * v32 * v34 * v60 * v61
          + v11 * v12 * v22 * v26 * v32 * v35 * v60 * v61
          + v11 * v13 * v22 * v26 * v32 * v35 * v60 * v61
          + v10 * v13 * v23 * v31 * v33 * v37 * v60 * v61
          + v11 * v13 * v23 * v31 * v33 * v37 * v60 * v61
          + v10 * v13 * v23 * v29 * v33 * v38 * v60 * v61
          + v11 * v13 * v23 * v29 * v33 * v38 * v60 * v61
          + v11 * v12 * v24 * v28 * v32 * v40 * v60 * v61
          + v11 * v13 * v24 * v28 * v32 * v40 * v60 * v61
          + v11 * v12 * v24 * v26 * v32 * v41 * v60 * v61
          + v11 * v13 * v24 * v26 * v32 * v41 * v60 * v61
          + v10 * v13 * v25 * v31 * v33 * v43 * v60 * v61
          + v11 * v13 * v25 * v31 * v33 * v43 * v60 * v61
          + v10 * v13 * v25 * v29 * v33 * v44 * v60 * v61
          + v11 * v13 * v25 * v29 * v33 * v44 * v60 * v61
          + v11 * v12 * v16 * v28 * v32 * v46 * v60 * v61
          + v11 * v13 * v16 * v28 * v32 * v46 * v60 * v61
          + v11 * v12 * v16 * v26 * v32 * v47 * v60 * v61
          + v11 * v13 * v16 * v26 * v32 * v47 * v60 * v61
          + v10 * v13 * v21 * v31 * v33 * v49 * v60 * v61
          + v11 * v13 * v21 * v31 * v33 * v49 * v60 * v61
          + v10 * v13 * v21 * v29 * v33 * v50 * v60 * v61
          + v11 * v13 * v21 * v29 * v33 * v50 * v60 * v61
          + v11 * v12 * v16 * v26 * v28 * v54 * v60 * v61
          + v11 * v13 * v16 * v26 * v28 * v54 * v60 * v61
          + v10 * v13 * v21 * v29 * v31 * v55 * v60 * v61
          + v11 * v13 * v21 * v29 * v31 * v55 * v60 * v61
          + v11 * v12 * v22 * v26 * v28 * v56 * v60 * v61
          + v11 * v13 * v22 * v26 * v28 * v56 * v60 * v61
          + v10 * v13 * v23 * v29 * v31 * v57 * v60 * v61
          + v11 * v13 * v23 * v29 * v31 * v57 * v60 * v61
          + v11 * v12 * v24 * v26 * v28 * v58 * v60 * v61
          + v11 * v13 * v24 * v26 * v28 * v58 * v60 * v61
          + v10 * v13 * v25 * v29 * v31 * v59 * v60 * v61
          + v11 * v13 * v25 * v29 * v31 * v59 * v60 * v61
          + sqrt(1 - v13 + v19) * v20 * v22 * v27 * v32 * v34 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v22 * v27 * v32 * v34 * v6 * v60 * v7
          + sqrt(1 - v13 + v19) * v20 * v22 * v26 * v32 * v36 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v22 * v26 * v32 * v36 * v6 * v60 * v7
          + sqrt(1 - v13 + v19) * v20 * v24 * v27 * v32 * v40 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v24 * v27 * v32 * v40 * v6 * v60 * v7
          + sqrt(1 - v13 + v19) * v20 * v24 * v26 * v32 * v42 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v24 * v26 * v32 * v42 * v6 * v60 * v7
          + v16 * sqrt(1 - v13 + v19) * v20 * v27 * v32 * v46 * v6 * v60 * v7
          + v16 * v19 * sqrt(1 - v12 + v20) * v27 * v32 * v46 * v6 * v60 * v7
          + v16 * sqrt(1 - v13 + v19) * v20 * v26 * v32 * v48 * v6 * v60 * v7
          + v16 * v19 * sqrt(1 - v12 + v20) * v26 * v32 * v48 * v6 * v60 * v7
          + v16 * sqrt(1 - v13 + v19) * v20 * v26 * v27 * v54 * v6 * v60 * v7
          + v16 * v19 * sqrt(1 - v12 + v20) * v26 * v27 * v54 * v6 * v60 * v7
          + sqrt(1 - v13 + v19) * v20 * v22 * v26 * v27 * v56 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v22 * v26 * v27 * v56 * v6 * v60 * v7
          + sqrt(1 - v13 + v19) * v20 * v24 * v26 * v27 * v58 * v6 * v60 * v7
          + v19 * sqrt(1 - v12 + v20) * v24 * v26 * v27 * v58 * v6 * v60 * v7
          + v12 * v22 * v27 * v32 * v34 * v6 * v61 * v7
          + v13 * v22 * v27 * v32 * v34 * v6 * v61 * v7
          + v12 * v22 * v26 * v32 * v36 * v6 * v61 * v7
          + v13 * v22 * v26 * v32 * v36 * v6 * v61 * v7
          + v12 * v24 * v27 * v32 * v40 * v6 * v61 * v7
          + v13 * v24 * v27 * v32 * v40 * v6 * v61 * v7
          + v12 * v24 * v26 * v32 * v42 * v6 * v61 * v7
          + v13 * v24 * v26 * v32 * v42 * v6 * v61 * v7
          + v12 * v16 * v27 * v32 * v46 * v6 * v61 * v7
          + v13 * v16 * v27 * v32 * v46 * v6 * v61 * v7
          + v12 * v16 * v26 * v32 * v48 * v6 * v61 * v7
          + v13 * v16 * v26 * v32 * v48 * v6 * v61 * v7
          + v12 * v16 * v26 * v27 * v54 * v6 * v61 * v7
          + v13 * v16 * v26 * v27 * v54 * v6 * v61 * v7
          + v12 * v22 * v26 * v27 * v56 * v6 * v61 * v7
          + v13 * v22 * v26 * v27 * v56 * v6 * v61 * v7
          + v12 * v24 * v26 * v27 * v58 * v6 * v61 * v7
          + v13 * v24 * v26 * v27 * v58 * v6 * v61 * v7
          + v13 * sqrt(1 - v11 + v17) * v18 * v23 * v33 * v37 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v23 * v33 * v37 * v8
          + v13 * sqrt(1 - v11 + v17) * v18 * v25 * v33 * v43 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v25 * v33 * v43 * v8
          + v13 * sqrt(1 - v11 + v17) * v18 * v21 * v33 * v49 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v21 * v33 * v49 * v8
          + v13 * sqrt(1 - v11 + v17) * v18 * v21 * v29 * v55 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v21 * v29 * v55 * v8
          + v13 * sqrt(1 - v11 + v17) * v18 * v23 * v29 * v57 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v23 * v29 * v57 * v8
          + v13 * sqrt(1 - v11 + v17) * v18 * v25 * v29 * v59 * v8
          + v13 * v17 * sqrt(1 - v10 + v18) * v25 * v29 * v59 * v8
          + v10 * v23 * v30 * v33 * v37 * v60 * v8 * v9
          + v11 * v23 * v30 * v33 * v37 * v60 * v8 * v9
          + v10 * v23 * v29 * v33 * v39 * v60 * v8 * v9
          + v11 * v23 * v29 * v33 * v39 * v60 * v8 * v9
          + v10 * v25 * v30 * v33 * v43 * v60 * v8 * v9
          + v11 * v25 * v30 * v33 * v43 * v60 * v8 * v9
          + v10 * v25 * v29 * v33 * v45 * v60 * v8 * v9
          + v11 * v25 * v29 * v33 * v45 * v60 * v8 * v9
          + v10 * v21 * v30 * v33 * v49 * v60 * v8 * v9
          + v11 * v21 * v30 * v33 * v49 * v60 * v8 * v9
          + v10 * v21 * v29 * v33 * v51 * v60 * v8 * v9
          + v11 * v21 * v29 * v33 * v51 * v60 * v8 * v9
          + v10 * v21 * v29 * v30 * v55 * v60 * v8 * v9
          + v11 * v21 * v29 * v30 * v55 * v60 * v8 * v9
          + v10 * v23 * v29 * v30 * v57 * v60 * v8 * v9
          + v11 * v23 * v29 * v30 * v57 * v60 * v8 * v9
          + v10 * v25 * v29 * v30 * v59 * v60 * v8 * v9
          + v11 * v25 * v29 * v30 * v59 * v60 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v23 * v30 * v33 * v37 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v23 * v30 * v33 * v37 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v23 * v29 * v33 * v39 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v23 * v29 * v33 * v39 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v25 * v30 * v33 * v43 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v25 * v30 * v33 * v43 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v25 * v29 * v33 * v45 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v25 * v29 * v33 * v45 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v21 * v30 * v33 * v49 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v21 * v30 * v33 * v49 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v21 * v29 * v33 * v51 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v21 * v29 * v33 * v51 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v21 * v29 * v30 * v55 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v21 * v29 * v30 * v55 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v23 * v29 * v30 * v57 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v23 * v29 * v30 * v57 * v61 * v8 * v9
          + sqrt(1 - v11 + v17) * v18 * v25 * v29 * v30 * v59 * v61 * v8 * v9
          + v17 * sqrt(1 - v10 + v18) * v25 * v29 * v30 * v59 * v61 * v8 * v9)
         / ((v10 + v11) * (v12 + v13) * sqrt(v52) * sqrt(v53));
}