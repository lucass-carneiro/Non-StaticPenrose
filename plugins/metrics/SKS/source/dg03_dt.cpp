#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::dllgSKS_03_dt(double t, double x, double y, double z) const noexcept
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
  const double v7{dsy_dt(1, t)};
  const double v8{dsx_dt(1, t)};
  const double v9{dsx_dt(2, t)};
  const double v10{d2sy_dt2(2, t)};
  const double v11{d2sx_dt2(2, t)};
  const double v12{d2sy_dt2(1, t)};
  const double v13{d2sx_dt2(1, t)};
  const double v14{pow(v6, 3)};
  const double v15{pow(v7, 2)};
  const double v16{pow(v8, 2)};
  const double v17{pow(v6, 2)};
  const double v18{pow(v9, 2)};
  const double v19{pow(v9, 3)};
  const double v20{pow(v7, 3)};
  const double v21{pow(v8, 3)};
  const double v22{dZ_dt(2, t, x, y, z)};
  const double v23{-v15};
  const double v24{-v16};
  const double v25{dZ_dt(1, t, x, y, z)};
  const double v26{-v17};
  const double v27{-v18};
  const double v28{dY_dt(2, t, x, y, z)};
  const double v29{dY_dt(1, t, x, y, z)};
  const double v30{dX_dt(2, t, x, y, z)};
  const double v31{dX_dt(1, t, x, y, z)};
  const double v32{l3_KS(a2, v2, v1, v0)};
  const double v33{l2_KS(a2, v2, v1, v0)};
  const double v34{l1_KS(a2, v2, v1, v0)};
  const double v35{l3_KS(a1, v5, v4, v3)};
  const double v36{l2_KS(a1, v5, v4, v3)};
  const double v37{l1_KS(a1, v5, v4, v3)};
  const double v38{H_KS(M2, a2, v2, v1, v0)};
  const double v39{H_KS(M1, a1, v5, v4, v3)};
  const double v40{dl3_KS_dZ(a2, v2, v1, v0)};
  const double v41{dl2_KS_dZ(a2, v2, v1, v0)};
  const double v42{dl1_KS_dZ(a2, v2, v1, v0)};
  const double v43{dl3_KS_dZ(a1, v5, v4, v3)};
  const double v44{dl2_KS_dZ(a1, v5, v4, v3)};
  const double v45{dl1_KS_dZ(a1, v5, v4, v3)};
  const double v46{dl3_KS_dY(a2, v2, v1, v0)};
  const double v47{dl2_KS_dY(a2, v2, v1, v0)};
  const double v48{dl1_KS_dY(a2, v2, v1, v0)};
  const double v49{dl3_KS_dY(a1, v5, v4, v3)};
  const double v50{dl2_KS_dY(a1, v5, v4, v3)};
  const double v51{dl1_KS_dY(a1, v5, v4, v3)};
  const double v52{dl3_KS_dX(a2, v2, v1, v0)};
  const double v53{dl2_KS_dX(a2, v2, v1, v0)};
  const double v54{dl1_KS_dX(a2, v2, v1, v0)};
  const double v55{dl3_KS_dX(a1, v5, v4, v3)};
  const double v56{dl2_KS_dX(a1, v5, v4, v3)};
  const double v57{dl1_KS_dX(a1, v5, v4, v3)};
  const double v58{1 + v23 + v24};
  const double v59{1 + v26 + v27};
  const double v60{dH_KS_dZ(M2, a2, v2, v1, v0)};
  const double v61{dH_KS_dZ(M1, a1, v5, v4, v3)};
  const double v62{dH_KS_dY(M2, a2, v2, v1, v0)};
  const double v63{dH_KS_dY(M1, a1, v5, v4, v3)};
  const double v64{dH_KS_dX(M2, a2, v2, v1, v0)};
  const double v65{dH_KS_dX(M1, a1, v5, v4, v3)};
  const double v66{sqrt(v58)};
  const double v67{sqrt(v59)};

  return (v10 * v18 * sqrt(1 - v16 + v23) * v24 * v32 * v33 * v38
          + v10 * v18 * v23 * sqrt(1 - v15 + v24) * v32 * v33 * v38
          + v11 * v17 * sqrt(1 - v16 + v23) * v24 * v32 * v34 * v38
          + v11 * v17 * v23 * sqrt(1 - v15 + v24) * v32 * v34 * v38
          + v22 * sqrt(1 - v16 + v23) * v24 * v38 * v40
          + v22 * v23 * sqrt(1 - v15 + v24) * v38 * v40
          + v14 * v22 * sqrt(1 - v16 + v23) * v24 * v33 * v38 * v40
          + v14 * v22 * v23 * sqrt(1 - v15 + v24) * v33 * v38 * v40
          + v19 * v22 * sqrt(1 - v16 + v23) * v24 * v34 * v38 * v40
          + v19 * v22 * v23 * sqrt(1 - v15 + v24) * v34 * v38 * v40
          + v14 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v41
          + v14 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v41
          + v19 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v42
          + v19 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v42
          + v25 * sqrt(1 - v18 + v26) * v27 * v39 * v43
          + v25 * v26 * sqrt(1 - v17 + v27) * v39 * v43
          + v20 * v25 * sqrt(1 - v18 + v26) * v27 * v36 * v39 * v43
          + v20 * v25 * v26 * sqrt(1 - v17 + v27) * v36 * v39 * v43
          + v21 * v25 * sqrt(1 - v18 + v26) * v27 * v37 * v39 * v43
          + v21 * v25 * v26 * sqrt(1 - v17 + v27) * v37 * v39 * v43
          + v20 * v25 * sqrt(1 - v18 + v26) * v27 * v35 * v39 * v44
          + v20 * v25 * v26 * sqrt(1 - v17 + v27) * v35 * v39 * v44
          + v21 * v25 * sqrt(1 - v18 + v26) * v27 * v35 * v39 * v45
          + v21 * v25 * v26 * sqrt(1 - v17 + v27) * v35 * v39 * v45
          + sqrt(1 - v16 + v23) * v24 * v28 * v38 * v46
          + v23 * sqrt(1 - v15 + v24) * v28 * v38 * v46
          + v14 * sqrt(1 - v16 + v23) * v24 * v28 * v33 * v38 * v46
          + v14 * v23 * sqrt(1 - v15 + v24) * v28 * v33 * v38 * v46
          + v19 * sqrt(1 - v16 + v23) * v24 * v28 * v34 * v38 * v46
          + v19 * v23 * sqrt(1 - v15 + v24) * v28 * v34 * v38 * v46
          + v14 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v38 * v47
          + v14 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v38 * v47
          + v19 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v38 * v48
          + v19 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v38 * v48
          + sqrt(1 - v18 + v26) * v27 * v29 * v39 * v49
          + v26 * sqrt(1 - v17 + v27) * v29 * v39 * v49
          + v20 * sqrt(1 - v18 + v26) * v27 * v29 * v36 * v39 * v49
          + v20 * v26 * sqrt(1 - v17 + v27) * v29 * v36 * v39 * v49
          + v21 * sqrt(1 - v18 + v26) * v27 * v29 * v37 * v39 * v49
          + v21 * v26 * sqrt(1 - v17 + v27) * v29 * v37 * v39 * v49
          + v20 * sqrt(1 - v18 + v26) * v27 * v29 * v35 * v39 * v50
          + v20 * v26 * sqrt(1 - v17 + v27) * v29 * v35 * v39 * v50
          + v21 * sqrt(1 - v18 + v26) * v27 * v29 * v35 * v39 * v51
          + v21 * v26 * sqrt(1 - v17 + v27) * v29 * v35 * v39 * v51
          + sqrt(1 - v16 + v23) * v24 * v30 * v38 * v52
          + v23 * sqrt(1 - v15 + v24) * v30 * v38 * v52
          + v14 * sqrt(1 - v16 + v23) * v24 * v30 * v33 * v38 * v52
          + v14 * v23 * sqrt(1 - v15 + v24) * v30 * v33 * v38 * v52
          + v19 * sqrt(1 - v16 + v23) * v24 * v30 * v34 * v38 * v52
          + v19 * v23 * sqrt(1 - v15 + v24) * v30 * v34 * v38 * v52
          + v14 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v38 * v53
          + v14 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v38 * v53
          + v19 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v38 * v54
          + v19 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v38 * v54
          + sqrt(1 - v18 + v26) * v27 * v31 * v39 * v55
          + v26 * sqrt(1 - v17 + v27) * v31 * v39 * v55
          + v20 * sqrt(1 - v18 + v26) * v27 * v31 * v36 * v39 * v55
          + v20 * v26 * sqrt(1 - v17 + v27) * v31 * v36 * v39 * v55
          + v21 * sqrt(1 - v18 + v26) * v27 * v31 * v37 * v39 * v55
          + v21 * v26 * sqrt(1 - v17 + v27) * v31 * v37 * v39 * v55
          + v20 * sqrt(1 - v18 + v26) * v27 * v31 * v35 * v39 * v56
          + v20 * v26 * sqrt(1 - v17 + v27) * v31 * v35 * v39 * v56
          + v21 * sqrt(1 - v18 + v26) * v27 * v31 * v35 * v39 * v57
          + v21 * v26 * sqrt(1 - v17 + v27) * v31 * v35 * v39 * v57
          + v10 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v6
          + v10 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v6
          + v18 * v22 * sqrt(1 - v16 + v23) * v24 * v33 * v38 * v40 * v6
          + v18 * v22 * v23 * sqrt(1 - v15 + v24) * v33 * v38 * v40 * v6
          + v18 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v41 * v6
          + v18 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v41 * v6
          + v18 * sqrt(1 - v16 + v23) * v24 * v28 * v33 * v38 * v46 * v6
          + v18 * v23 * sqrt(1 - v15 + v24) * v28 * v33 * v38 * v46 * v6
          + v18 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v38 * v47 * v6
          + v18 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v38 * v47 * v6
          + v18 * sqrt(1 - v16 + v23) * v24 * v30 * v33 * v38 * v52 * v6
          + v18 * v23 * sqrt(1 - v15 + v24) * v30 * v33 * v38 * v52 * v6
          + v18 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v38 * v53 * v6
          + v18 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v38 * v53 * v6
          + v22 * sqrt(1 - v16 + v23) * v24 * v32 * v60
          + v22 * v23 * sqrt(1 - v15 + v24) * v32 * v60
          + v14 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v33 * v60
          + v14 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v33 * v60
          + v19 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v34 * v60
          + v19 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v34 * v60
          + v18 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v33 * v6 * v60
          + v18 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v33 * v6 * v60
          + v25 * sqrt(1 - v18 + v26) * v27 * v35 * v61
          + v25 * v26 * sqrt(1 - v17 + v27) * v35 * v61
          + v20 * v25 * sqrt(1 - v18 + v26) * v27 * v35 * v36 * v61
          + v20 * v25 * v26 * sqrt(1 - v17 + v27) * v35 * v36 * v61
          + v21 * v25 * sqrt(1 - v18 + v26) * v27 * v35 * v37 * v61
          + v21 * v25 * v26 * sqrt(1 - v17 + v27) * v35 * v37 * v61
          + sqrt(1 - v16 + v23) * v24 * v28 * v32 * v62
          + v23 * sqrt(1 - v15 + v24) * v28 * v32 * v62
          + v14 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v33 * v62
          + v14 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v33 * v62
          + v19 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v34 * v62
          + v19 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v34 * v62
          + v18 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v33 * v6 * v62
          + v18 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v33 * v6 * v62
          + sqrt(1 - v18 + v26) * v27 * v29 * v35 * v63
          + v26 * sqrt(1 - v17 + v27) * v29 * v35 * v63
          + v20 * sqrt(1 - v18 + v26) * v27 * v29 * v35 * v36 * v63
          + v20 * v26 * sqrt(1 - v17 + v27) * v29 * v35 * v36 * v63
          + v21 * sqrt(1 - v18 + v26) * v27 * v29 * v35 * v37 * v63
          + v21 * v26 * sqrt(1 - v17 + v27) * v29 * v35 * v37 * v63
          + sqrt(1 - v16 + v23) * v24 * v30 * v32 * v64
          + v23 * sqrt(1 - v15 + v24) * v30 * v32 * v64
          + v14 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v33 * v64
          + v14 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v33 * v64
          + v19 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v34 * v64
          + v19 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v34 * v64
          + v18 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v33 * v6 * v64
          + v18 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v33 * v6 * v64
          + sqrt(1 - v18 + v26) * v27 * v31 * v35 * v65
          + v26 * sqrt(1 - v17 + v27) * v31 * v35 * v65
          + v20 * sqrt(1 - v18 + v26) * v27 * v31 * v35 * v36 * v65
          + v20 * v26 * sqrt(1 - v17 + v27) * v31 * v35 * v36 * v65
          + v21 * sqrt(1 - v18 + v26) * v27 * v31 * v35 * v37 * v65
          + v21 * v26 * sqrt(1 - v17 + v27) * v31 * v35 * v37 * v65 - v10 * v32 * v33 * v38 * v66
          + v10 * v15 * v32 * v33 * v38 * v66 + v10 * v16 * v32 * v33 * v38 * v66
          + v10 * v18 * v32 * v33 * v38 * v66 - v11 * v32 * v34 * v38 * v66
          + v11 * v15 * v32 * v34 * v38 * v66 + v11 * v16 * v32 * v34 * v38 * v66
          + v11 * v17 * v32 * v34 * v38 * v66 + v22 * v38 * v40 * v66
          + v15 * v17 * v22 * v38 * v40 * v66 + v16 * v17 * v22 * v38 * v40 * v66
          + v15 * v18 * v22 * v38 * v40 * v66 + v16 * v18 * v22 * v38 * v40 * v66
          + v22 * v26 * v38 * v40 * v66 + v22 * v27 * v38 * v40 * v66
          + v14 * v22 * v33 * v38 * v40 * v66 + v19 * v22 * v34 * v38 * v40 * v66
          + v14 * v22 * v32 * v38 * v41 * v66 + v19 * v22 * v32 * v38 * v42 * v66
          + v28 * v38 * v46 * v66 + v15 * v17 * v28 * v38 * v46 * v66
          + v16 * v17 * v28 * v38 * v46 * v66 + v15 * v18 * v28 * v38 * v46 * v66
          + v16 * v18 * v28 * v38 * v46 * v66 + v26 * v28 * v38 * v46 * v66
          + v27 * v28 * v38 * v46 * v66 + v14 * v28 * v33 * v38 * v46 * v66
          + v19 * v28 * v34 * v38 * v46 * v66 + v14 * v28 * v32 * v38 * v47 * v66
          + v19 * v28 * v32 * v38 * v48 * v66 + v30 * v38 * v52 * v66
          + v15 * v17 * v30 * v38 * v52 * v66 + v16 * v17 * v30 * v38 * v52 * v66
          + v15 * v18 * v30 * v38 * v52 * v66 + v16 * v18 * v30 * v38 * v52 * v66
          + v26 * v30 * v38 * v52 * v66 + v27 * v30 * v38 * v52 * v66
          + v14 * v30 * v33 * v38 * v52 * v66 + v19 * v30 * v34 * v38 * v52 * v66
          + v14 * v30 * v32 * v38 * v53 * v66 + v19 * v30 * v32 * v38 * v54 * v66
          + v10 * v32 * v38 * v6 * v66 - v22 * v33 * v38 * v40 * v6 * v66
          + v15 * v22 * v33 * v38 * v40 * v6 * v66 + v16 * v22 * v33 * v38 * v40 * v6 * v66
          + v18 * v22 * v33 * v38 * v40 * v6 * v66 - v22 * v32 * v38 * v41 * v6 * v66
          + v15 * v22 * v32 * v38 * v41 * v6 * v66 + v16 * v22 * v32 * v38 * v41 * v6 * v66
          + v18 * v22 * v32 * v38 * v41 * v6 * v66 - v28 * v33 * v38 * v46 * v6 * v66
          + v15 * v28 * v33 * v38 * v46 * v6 * v66 + v16 * v28 * v33 * v38 * v46 * v6 * v66
          + v18 * v28 * v33 * v38 * v46 * v6 * v66 - v28 * v32 * v38 * v47 * v6 * v66
          + v15 * v28 * v32 * v38 * v47 * v6 * v66 + v16 * v28 * v32 * v38 * v47 * v6 * v66
          + v18 * v28 * v32 * v38 * v47 * v6 * v66 - v30 * v33 * v38 * v52 * v6 * v66
          + v15 * v30 * v33 * v38 * v52 * v6 * v66 + v16 * v30 * v33 * v38 * v52 * v6 * v66
          + v18 * v30 * v33 * v38 * v52 * v6 * v66 - v30 * v32 * v38 * v53 * v6 * v66
          + v15 * v30 * v32 * v38 * v53 * v6 * v66 + v16 * v30 * v32 * v38 * v53 * v6 * v66
          + v18 * v30 * v32 * v38 * v53 * v6 * v66 + v22 * v32 * v60 * v66
          + v15 * v17 * v22 * v32 * v60 * v66 + v16 * v17 * v22 * v32 * v60 * v66
          + v15 * v18 * v22 * v32 * v60 * v66 + v16 * v18 * v22 * v32 * v60 * v66
          + v22 * v26 * v32 * v60 * v66 + v22 * v27 * v32 * v60 * v66
          + v14 * v22 * v32 * v33 * v60 * v66 + v19 * v22 * v32 * v34 * v60 * v66
          - v22 * v32 * v33 * v6 * v60 * v66 + v15 * v22 * v32 * v33 * v6 * v60 * v66
          + v16 * v22 * v32 * v33 * v6 * v60 * v66 + v18 * v22 * v32 * v33 * v6 * v60 * v66
          + v28 * v32 * v62 * v66 + v15 * v17 * v28 * v32 * v62 * v66
          + v16 * v17 * v28 * v32 * v62 * v66 + v15 * v18 * v28 * v32 * v62 * v66
          + v16 * v18 * v28 * v32 * v62 * v66 + v26 * v28 * v32 * v62 * v66
          + v27 * v28 * v32 * v62 * v66 + v14 * v28 * v32 * v33 * v62 * v66
          + v19 * v28 * v32 * v34 * v62 * v66 - v28 * v32 * v33 * v6 * v62 * v66
          + v15 * v28 * v32 * v33 * v6 * v62 * v66 + v16 * v28 * v32 * v33 * v6 * v62 * v66
          + v18 * v28 * v32 * v33 * v6 * v62 * v66 + v30 * v32 * v64 * v66
          + v15 * v17 * v30 * v32 * v64 * v66 + v16 * v17 * v30 * v32 * v64 * v66
          + v15 * v18 * v30 * v32 * v64 * v66 + v16 * v18 * v30 * v32 * v64 * v66
          + v26 * v30 * v32 * v64 * v66 + v27 * v30 * v32 * v64 * v66
          + v14 * v30 * v32 * v33 * v64 * v66 + v19 * v30 * v32 * v34 * v64 * v66
          - v30 * v32 * v33 * v6 * v64 * v66 + v15 * v30 * v32 * v33 * v6 * v64 * v66
          + v16 * v30 * v32 * v33 * v6 * v64 * v66 + v18 * v30 * v32 * v33 * v6 * v64 * v66
          - v12 * v35 * v36 * v39 * v67 + v12 * v16 * v35 * v36 * v39 * v67
          + v12 * v17 * v35 * v36 * v39 * v67 + v12 * v18 * v35 * v36 * v39 * v67
          + v12 * v17 * v24 * v35 * v36 * v39 * v67 + v12 * v18 * v24 * v35 * v36 * v39 * v67
          - v13 * v35 * v37 * v39 * v67 + v13 * v15 * v35 * v37 * v39 * v67
          + v13 * v17 * v35 * v37 * v39 * v67 + v13 * v18 * v35 * v37 * v39 * v67
          + v13 * v17 * v23 * v35 * v37 * v39 * v67 + v13 * v18 * v23 * v35 * v37 * v39 * v67
          + v25 * v39 * v43 * v67 + v15 * v17 * v25 * v39 * v43 * v67
          + v16 * v17 * v25 * v39 * v43 * v67 + v15 * v18 * v25 * v39 * v43 * v67
          + v16 * v18 * v25 * v39 * v43 * v67 + v23 * v25 * v39 * v43 * v67
          + v24 * v25 * v39 * v43 * v67 + v20 * v25 * v36 * v39 * v43 * v67
          + v21 * v25 * v37 * v39 * v43 * v67 + v20 * v25 * v35 * v39 * v44 * v67
          + v21 * v25 * v35 * v39 * v45 * v67 + v29 * v39 * v49 * v67
          + v15 * v17 * v29 * v39 * v49 * v67 + v16 * v17 * v29 * v39 * v49 * v67
          + v15 * v18 * v29 * v39 * v49 * v67 + v16 * v18 * v29 * v39 * v49 * v67
          + v23 * v29 * v39 * v49 * v67 + v24 * v29 * v39 * v49 * v67
          + v20 * v29 * v36 * v39 * v49 * v67 + v21 * v29 * v37 * v39 * v49 * v67
          + v20 * v29 * v35 * v39 * v50 * v67 + v21 * v29 * v35 * v39 * v51 * v67
          + v31 * v39 * v55 * v67 + v15 * v17 * v31 * v39 * v55 * v67
          + v16 * v17 * v31 * v39 * v55 * v67 + v15 * v18 * v31 * v39 * v55 * v67
          + v16 * v18 * v31 * v39 * v55 * v67 + v23 * v31 * v39 * v55 * v67
          + v24 * v31 * v39 * v55 * v67 + v20 * v31 * v36 * v39 * v55 * v67
          + v21 * v31 * v37 * v39 * v55 * v67 + v20 * v31 * v35 * v39 * v56 * v67
          + v21 * v31 * v35 * v39 * v57 * v67 + v25 * v35 * v61 * v67
          + v15 * v17 * v25 * v35 * v61 * v67 + v16 * v17 * v25 * v35 * v61 * v67
          + v15 * v18 * v25 * v35 * v61 * v67 + v16 * v18 * v25 * v35 * v61 * v67
          + v23 * v25 * v35 * v61 * v67 + v24 * v25 * v35 * v61 * v67
          + v20 * v25 * v35 * v36 * v61 * v67 + v21 * v25 * v35 * v37 * v61 * v67
          + v29 * v35 * v63 * v67 + v15 * v17 * v29 * v35 * v63 * v67
          + v16 * v17 * v29 * v35 * v63 * v67 + v15 * v18 * v29 * v35 * v63 * v67
          + v16 * v18 * v29 * v35 * v63 * v67 + v23 * v29 * v35 * v63 * v67
          + v24 * v29 * v35 * v63 * v67 + v20 * v29 * v35 * v36 * v63 * v67
          + v21 * v29 * v35 * v37 * v63 * v67 + v31 * v35 * v65 * v67
          + v15 * v17 * v31 * v35 * v65 * v67 + v16 * v17 * v31 * v35 * v65 * v67
          + v15 * v18 * v31 * v35 * v65 * v67 + v16 * v18 * v31 * v35 * v65 * v67
          + v23 * v31 * v35 * v65 * v67 + v24 * v31 * v35 * v65 * v67
          + v20 * v31 * v35 * v36 * v65 * v67 + v21 * v31 * v35 * v37 * v65 * v67
          + v12 * sqrt(1 - v18 + v26) * v27 * v35 * v39 * v7
          + v12 * v26 * sqrt(1 - v17 + v27) * v35 * v39 * v7 + v12 * v35 * v39 * v67 * v7
          - v25 * v36 * v39 * v43 * v67 * v7 + v16 * v25 * v36 * v39 * v43 * v67 * v7
          + v17 * v25 * v36 * v39 * v43 * v67 * v7 + v18 * v25 * v36 * v39 * v43 * v67 * v7
          + v17 * v24 * v25 * v36 * v39 * v43 * v67 * v7
          + v18 * v24 * v25 * v36 * v39 * v43 * v67 * v7 - v25 * v35 * v39 * v44 * v67 * v7
          + v16 * v25 * v35 * v39 * v44 * v67 * v7 + v17 * v25 * v35 * v39 * v44 * v67 * v7
          + v18 * v25 * v35 * v39 * v44 * v67 * v7 + v17 * v24 * v25 * v35 * v39 * v44 * v67 * v7
          + v18 * v24 * v25 * v35 * v39 * v44 * v67 * v7 - v29 * v36 * v39 * v49 * v67 * v7
          + v16 * v29 * v36 * v39 * v49 * v67 * v7 + v17 * v29 * v36 * v39 * v49 * v67 * v7
          + v18 * v29 * v36 * v39 * v49 * v67 * v7 + v17 * v24 * v29 * v36 * v39 * v49 * v67 * v7
          + v18 * v24 * v29 * v36 * v39 * v49 * v67 * v7 - v29 * v35 * v39 * v50 * v67 * v7
          + v16 * v29 * v35 * v39 * v50 * v67 * v7 + v17 * v29 * v35 * v39 * v50 * v67 * v7
          + v18 * v29 * v35 * v39 * v50 * v67 * v7 + v17 * v24 * v29 * v35 * v39 * v50 * v67 * v7
          + v18 * v24 * v29 * v35 * v39 * v50 * v67 * v7 - v31 * v36 * v39 * v55 * v67 * v7
          + v16 * v31 * v36 * v39 * v55 * v67 * v7 + v17 * v31 * v36 * v39 * v55 * v67 * v7
          + v18 * v31 * v36 * v39 * v55 * v67 * v7 + v17 * v24 * v31 * v36 * v39 * v55 * v67 * v7
          + v18 * v24 * v31 * v36 * v39 * v55 * v67 * v7 - v31 * v35 * v39 * v56 * v67 * v7
          + v16 * v31 * v35 * v39 * v56 * v67 * v7 + v17 * v31 * v35 * v39 * v56 * v67 * v7
          + v18 * v31 * v35 * v39 * v56 * v67 * v7 + v17 * v24 * v31 * v35 * v39 * v56 * v67 * v7
          + v18 * v24 * v31 * v35 * v39 * v56 * v67 * v7 - v25 * v35 * v36 * v61 * v67 * v7
          + v16 * v25 * v35 * v36 * v61 * v67 * v7 + v17 * v25 * v35 * v36 * v61 * v67 * v7
          + v18 * v25 * v35 * v36 * v61 * v67 * v7 + v17 * v24 * v25 * v35 * v36 * v61 * v67 * v7
          + v18 * v24 * v25 * v35 * v36 * v61 * v67 * v7 - v29 * v35 * v36 * v63 * v67 * v7
          + v16 * v29 * v35 * v36 * v63 * v67 * v7 + v17 * v29 * v35 * v36 * v63 * v67 * v7
          + v18 * v29 * v35 * v36 * v63 * v67 * v7 + v17 * v24 * v29 * v35 * v36 * v63 * v67 * v7
          + v18 * v24 * v29 * v35 * v36 * v63 * v67 * v7 - v31 * v35 * v36 * v65 * v67 * v7
          + v16 * v31 * v35 * v36 * v65 * v67 * v7 + v17 * v31 * v35 * v36 * v65 * v67 * v7
          + v18 * v31 * v35 * v36 * v65 * v67 * v7 + v17 * v24 * v31 * v35 * v36 * v65 * v67 * v7
          + v18 * v24 * v31 * v35 * v36 * v65 * v67 * v7
          + v13 * sqrt(1 - v18 + v26) * v27 * v35 * v39 * v8
          + v13 * v26 * sqrt(1 - v17 + v27) * v35 * v39 * v8 + v13 * v35 * v39 * v67 * v8
          - v25 * v37 * v39 * v43 * v67 * v8 + v15 * v25 * v37 * v39 * v43 * v67 * v8
          + v17 * v25 * v37 * v39 * v43 * v67 * v8 + v18 * v25 * v37 * v39 * v43 * v67 * v8
          + v17 * v23 * v25 * v37 * v39 * v43 * v67 * v8
          + v18 * v23 * v25 * v37 * v39 * v43 * v67 * v8 - v25 * v35 * v39 * v45 * v67 * v8
          + v15 * v25 * v35 * v39 * v45 * v67 * v8 + v17 * v25 * v35 * v39 * v45 * v67 * v8
          + v18 * v25 * v35 * v39 * v45 * v67 * v8 + v17 * v23 * v25 * v35 * v39 * v45 * v67 * v8
          + v18 * v23 * v25 * v35 * v39 * v45 * v67 * v8 - v29 * v37 * v39 * v49 * v67 * v8
          + v15 * v29 * v37 * v39 * v49 * v67 * v8 + v17 * v29 * v37 * v39 * v49 * v67 * v8
          + v18 * v29 * v37 * v39 * v49 * v67 * v8 + v17 * v23 * v29 * v37 * v39 * v49 * v67 * v8
          + v18 * v23 * v29 * v37 * v39 * v49 * v67 * v8 - v29 * v35 * v39 * v51 * v67 * v8
          + v15 * v29 * v35 * v39 * v51 * v67 * v8 + v17 * v29 * v35 * v39 * v51 * v67 * v8
          + v18 * v29 * v35 * v39 * v51 * v67 * v8 + v17 * v23 * v29 * v35 * v39 * v51 * v67 * v8
          + v18 * v23 * v29 * v35 * v39 * v51 * v67 * v8 - v31 * v37 * v39 * v55 * v67 * v8
          + v15 * v31 * v37 * v39 * v55 * v67 * v8 + v17 * v31 * v37 * v39 * v55 * v67 * v8
          + v18 * v31 * v37 * v39 * v55 * v67 * v8 + v17 * v23 * v31 * v37 * v39 * v55 * v67 * v8
          + v18 * v23 * v31 * v37 * v39 * v55 * v67 * v8 - v31 * v35 * v39 * v57 * v67 * v8
          + v15 * v31 * v35 * v39 * v57 * v67 * v8 + v17 * v31 * v35 * v39 * v57 * v67 * v8
          + v18 * v31 * v35 * v39 * v57 * v67 * v8 + v17 * v23 * v31 * v35 * v39 * v57 * v67 * v8
          + v18 * v23 * v31 * v35 * v39 * v57 * v67 * v8 - v25 * v35 * v37 * v61 * v67 * v8
          + v15 * v25 * v35 * v37 * v61 * v67 * v8 + v17 * v25 * v35 * v37 * v61 * v67 * v8
          + v18 * v25 * v35 * v37 * v61 * v67 * v8 + v17 * v23 * v25 * v35 * v37 * v61 * v67 * v8
          + v18 * v23 * v25 * v35 * v37 * v61 * v67 * v8 - v29 * v35 * v37 * v63 * v67 * v8
          + v15 * v29 * v35 * v37 * v63 * v67 * v8 + v17 * v29 * v35 * v37 * v63 * v67 * v8
          + v18 * v29 * v35 * v37 * v63 * v67 * v8 + v17 * v23 * v29 * v35 * v37 * v63 * v67 * v8
          + v18 * v23 * v29 * v35 * v37 * v63 * v67 * v8 - v31 * v35 * v37 * v65 * v67 * v8
          + v15 * v31 * v35 * v37 * v65 * v67 * v8 + v17 * v31 * v35 * v37 * v65 * v67 * v8
          + v18 * v31 * v35 * v37 * v65 * v67 * v8 + v17 * v23 * v31 * v35 * v37 * v65 * v67 * v8
          + v18 * v23 * v31 * v35 * v37 * v65 * v67 * v8 - v13 * v35 * v36 * v39 * v67 * v7 * v8
          + v13 * v17 * v35 * v36 * v39 * v67 * v7 * v8
          + v13 * v18 * v35 * v36 * v39 * v67 * v7 * v8 - v12 * v35 * v37 * v39 * v67 * v7 * v8
          + v12 * v17 * v35 * v37 * v39 * v67 * v7 * v8
          + v12 * v18 * v35 * v37 * v39 * v67 * v7 * v8
          + v11 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v9
          + v11 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v9
          + v17 * v22 * sqrt(1 - v16 + v23) * v24 * v34 * v38 * v40 * v9
          + v17 * v22 * v23 * sqrt(1 - v15 + v24) * v34 * v38 * v40 * v9
          + v17 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v38 * v42 * v9
          + v17 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v38 * v42 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v28 * v34 * v38 * v46 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v28 * v34 * v38 * v46 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v38 * v48 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v38 * v48 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v30 * v34 * v38 * v52 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v30 * v34 * v38 * v52 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v38 * v54 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v38 * v54 * v9
          + v17 * v22 * sqrt(1 - v16 + v23) * v24 * v32 * v34 * v60 * v9
          + v17 * v22 * v23 * sqrt(1 - v15 + v24) * v32 * v34 * v60 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v28 * v32 * v34 * v62 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v28 * v32 * v34 * v62 * v9
          + v17 * sqrt(1 - v16 + v23) * v24 * v30 * v32 * v34 * v64 * v9
          + v17 * v23 * sqrt(1 - v15 + v24) * v30 * v32 * v34 * v64 * v9
          + v11 * v32 * v38 * v66 * v9 - v22 * v34 * v38 * v40 * v66 * v9
          + v15 * v22 * v34 * v38 * v40 * v66 * v9 + v16 * v22 * v34 * v38 * v40 * v66 * v9
          + v17 * v22 * v34 * v38 * v40 * v66 * v9 - v22 * v32 * v38 * v42 * v66 * v9
          + v15 * v22 * v32 * v38 * v42 * v66 * v9 + v16 * v22 * v32 * v38 * v42 * v66 * v9
          + v17 * v22 * v32 * v38 * v42 * v66 * v9 - v28 * v34 * v38 * v46 * v66 * v9
          + v15 * v28 * v34 * v38 * v46 * v66 * v9 + v16 * v28 * v34 * v38 * v46 * v66 * v9
          + v17 * v28 * v34 * v38 * v46 * v66 * v9 - v28 * v32 * v38 * v48 * v66 * v9
          + v15 * v28 * v32 * v38 * v48 * v66 * v9 + v16 * v28 * v32 * v38 * v48 * v66 * v9
          + v17 * v28 * v32 * v38 * v48 * v66 * v9 - v30 * v34 * v38 * v52 * v66 * v9
          + v15 * v30 * v34 * v38 * v52 * v66 * v9 + v16 * v30 * v34 * v38 * v52 * v66 * v9
          + v17 * v30 * v34 * v38 * v52 * v66 * v9 - v30 * v32 * v38 * v54 * v66 * v9
          + v15 * v30 * v32 * v38 * v54 * v66 * v9 + v16 * v30 * v32 * v38 * v54 * v66 * v9
          + v17 * v30 * v32 * v38 * v54 * v66 * v9 - v11 * v32 * v33 * v38 * v6 * v66 * v9
          + v11 * v15 * v32 * v33 * v38 * v6 * v66 * v9
          + v11 * v16 * v32 * v33 * v38 * v6 * v66 * v9 - v10 * v32 * v34 * v38 * v6 * v66 * v9
          + v10 * v15 * v32 * v34 * v38 * v6 * v66 * v9
          + v10 * v16 * v32 * v34 * v38 * v6 * v66 * v9 - v22 * v32 * v34 * v60 * v66 * v9
          + v15 * v22 * v32 * v34 * v60 * v66 * v9 + v16 * v22 * v32 * v34 * v60 * v66 * v9
          + v17 * v22 * v32 * v34 * v60 * v66 * v9 - v28 * v32 * v34 * v62 * v66 * v9
          + v15 * v28 * v32 * v34 * v62 * v66 * v9 + v16 * v28 * v32 * v34 * v62 * v66 * v9
          + v17 * v28 * v32 * v34 * v62 * v66 * v9 - v30 * v32 * v34 * v64 * v66 * v9
          + v15 * v30 * v32 * v34 * v64 * v66 * v9 + v16 * v30 * v32 * v34 * v64 * v66 * v9
          + v17 * v30 * v32 * v34 * v64 * v66 * v9)
         / ((-1 + v15 + v16) * (-1 + v17 + v18) * sqrt(v58) * sqrt(v59));
}