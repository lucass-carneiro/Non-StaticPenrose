#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_22(double t, double x, double y, double z) const noexcept -> double {

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
  const double v12{pow(v6, 4)};
  const double v13{pow(v8, 6)};
  const double v14{pow(v8, 4)};
  const double v15{pow(v9, 2)};
  const double v16{pow(v8, 2)};
  const double v17{pow(v9, 4)};
  const double v18{pow(v9, 6)};
  const double v19{pow(v6, 3)};
  const double v20{pow(v7, 3)};
  const double v21{pow(v7, 4)};
  const double v22{pow(v6, 6)};
  const double v23{pow(v8, 3)};
  const double v24{pow(v9, 3)};
  const double v25{pow(v8, 5)};
  const double v26{pow(v9, 5)};
  const double v27{pow(v6, 5)};
  const double v28{pow(v7, 5)};
  const double v29{pow(v7, 6)};
  const double v30{-v10};
  const double v31{-v11};
  const double v32{-v16};
  const double v33{-v15};
  const double v34{l1_KS(a2, v2, v1, v0)};
  const double v35{l2_KS(a2, v2, v1, v0)};
  const double v36{l1_KS(a1, v5, v4, v3)};
  const double v37{l2_KS(a1, v5, v4, v3)};
  const double v38{H_KS(M2, a2, v2, v1, v0)};
  const double v39{H_KS(M1, a1, v5, v4, v3)};
  const double v40{pow(v35, 2)};
  const double v41{pow(v34, 2)};
  const double v42{pow(v37, 2)};
  const double v43{pow(v36, 2)};
  const double v44{sqrt(1 + v30 + v31)};
  const double v45{sqrt(1 + v32 + v33)};

  return (-2 * v10 * v11 * v13 - v12 * v13 + 3 * v11 * v12 * v13 + 2 * v10 * v11 * v14 + v12 * v14
          - 3 * v11 * v12 * v14 - 6 * v10 * v11 * v14 * v15 - 3 * v12 * v14 * v15
          + 9 * v11 * v12 * v14 * v15 + 4 * v10 * v11 * v15 * v16 + 2 * v12 * v15 * v16
          - 6 * v11 * v12 * v15 * v16 + 2 * v10 * v11 * v17 + v12 * v17 - 3 * v11 * v12 * v17
          - 6 * v10 * v11 * v16 * v17 - 3 * v12 * v16 * v17 + 9 * v11 * v12 * v16 * v17
          - 2 * v10 * v11 * v18 - v12 * v18 + 3 * v11 * v12 * v18 - v13 * v21 + 3 * v10 * v13 * v21
          + v14 * v21 - 3 * v10 * v14 * v21 - 3 * v14 * v15 * v21 + 9 * v10 * v14 * v15 * v21
          + 2 * v15 * v16 * v21 - 6 * v10 * v15 * v16 * v21 + v17 * v21 - 3 * v10 * v17 * v21
          - 3 * v16 * v17 * v21 + 9 * v10 * v16 * v17 * v21 - v18 * v21 + 3 * v10 * v18 * v21
          + v13 * v22 - v14 * v22 + 3 * v14 * v15 * v22 - 2 * v15 * v16 * v22 - v17 * v22
          + 3 * v16 * v17 * v22 + v18 * v22 + v13 * v29 - v14 * v29 + 3 * v14 * v15 * v29
          - 2 * v15 * v16 * v29 - v17 * v29 + 3 * v16 * v17 * v29 + v18 * v29
          - 2 * v11 * v12 * v13 * v38 + 2 * v11 * v12 * v14 * v38 - 6 * v11 * v12 * v14 * v15 * v38
          + 4 * v11 * v12 * v15 * v16 * v38 + 2 * v11 * v12 * v17 * v38
          - 6 * v11 * v12 * v16 * v17 * v38 - 2 * v11 * v12 * v18 * v38 + v10 * v14 * v21 * v38
          - 3 * v10 * v14 * v15 * v21 * v38 + 2 * v10 * v15 * v16 * v21 * v38
          + v10 * v17 * v21 * v38 - 3 * v10 * v16 * v17 * v21 * v38 - v13 * v22 * v38
          + v14 * v22 * v38 - 3 * v14 * v15 * v22 * v38 + 2 * v15 * v16 * v22 * v38
          + v17 * v22 * v38 - 3 * v16 * v17 * v22 * v38 - v18 * v22 * v38 + v13 * v21 * v30 * v38
          + v18 * v21 * v30 * v38 + 2 * v10 * v13 * v20 * v34 * v38
          - 2 * v10 * v14 * v20 * v34 * v38 + 6 * v10 * v14 * v15 * v20 * v34 * v38
          - 4 * v10 * v15 * v16 * v20 * v34 * v38 - 2 * v10 * v17 * v20 * v34 * v38
          + 6 * v10 * v16 * v17 * v20 * v34 * v38 + 2 * v10 * v18 * v20 * v34 * v38
          + 2 * v11 * v13 * v19 * v35 * v38 - 2 * v11 * v14 * v19 * v35 * v38
          + 6 * v11 * v14 * v15 * v19 * v35 * v38 - 4 * v11 * v15 * v16 * v19 * v35 * v38
          - 2 * v11 * v17 * v19 * v35 * v38 + 6 * v11 * v16 * v17 * v19 * v35 * v38
          + 2 * v11 * v18 * v19 * v35 * v38 + 2 * v13 * v27 * v35 * v38 - 2 * v14 * v27 * v35 * v38
          + 6 * v14 * v15 * v27 * v35 * v38 - 4 * v15 * v16 * v27 * v35 * v38
          - 2 * v17 * v27 * v35 * v38 + 6 * v16 * v17 * v27 * v35 * v38 + 2 * v18 * v27 * v35 * v38
          - 2 * v13 * v19 * v20 * v34 * v35 * v38 + 2 * v14 * v19 * v20 * v34 * v35 * v38
          - 6 * v14 * v15 * v19 * v20 * v34 * v35 * v38
          + 4 * v15 * v16 * v19 * v20 * v34 * v35 * v38 + 2 * v17 * v19 * v20 * v34 * v35 * v38
          - 6 * v16 * v17 * v19 * v20 * v34 * v35 * v38 - 2 * v18 * v19 * v20 * v34 * v35 * v38
          + 2 * v10 * v11 * v13 * v39 + v12 * v13 * v39 - 3 * v11 * v12 * v13 * v39
          + 4 * v10 * v11 * v14 * v15 * v39 + 2 * v12 * v14 * v15 * v39
          - 6 * v11 * v12 * v14 * v15 * v39 + 2 * v10 * v11 * v16 * v17 * v39
          + v12 * v16 * v17 * v39 - 3 * v11 * v12 * v16 * v17 * v39 + v13 * v21 * v39
          - 3 * v10 * v13 * v21 * v39 + 2 * v14 * v15 * v21 * v39 - 6 * v10 * v14 * v15 * v21 * v39
          + v16 * v17 * v21 * v39 - 3 * v10 * v16 * v17 * v21 * v39 - v13 * v22 * v39
          - 2 * v14 * v15 * v22 * v39 - v13 * v29 * v39 - 2 * v14 * v15 * v29 * v39
          + v17 * v22 * v32 * v39 + v17 * v29 * v32 * v39 - 4 * v10 * v11 * v16 * v24 * v36 * v39
          - 2 * v12 * v16 * v24 * v36 * v39 + 6 * v11 * v12 * v16 * v24 * v36 * v39
          - 2 * v16 * v21 * v24 * v36 * v39 + 6 * v10 * v16 * v21 * v24 * v36 * v39
          + 2 * v16 * v22 * v24 * v36 * v39 + 2 * v16 * v24 * v29 * v36 * v39
          - 4 * v10 * v11 * v15 * v23 * v37 * v39 - 2 * v12 * v15 * v23 * v37 * v39
          + 6 * v11 * v12 * v15 * v23 * v37 * v39 - 2 * v15 * v21 * v23 * v37 * v39
          + 6 * v10 * v15 * v21 * v23 * v37 * v39 + 2 * v15 * v22 * v23 * v37 * v39
          - 4 * v10 * v11 * v25 * v37 * v39 - 2 * v12 * v25 * v37 * v39
          + 6 * v11 * v12 * v25 * v37 * v39 - 2 * v21 * v25 * v37 * v39
          + 6 * v10 * v21 * v25 * v37 * v39 + 2 * v22 * v25 * v37 * v39
          + 2 * v15 * v23 * v29 * v37 * v39 + 2 * v25 * v29 * v37 * v39
          + 4 * v10 * v11 * v23 * v24 * v36 * v37 * v39 + 2 * v12 * v23 * v24 * v36 * v37 * v39
          - 6 * v11 * v12 * v23 * v24 * v36 * v37 * v39 + 2 * v21 * v23 * v24 * v36 * v37 * v39
          - 6 * v10 * v21 * v23 * v24 * v36 * v37 * v39 - 2 * v22 * v23 * v24 * v36 * v37 * v39
          - 2 * v23 * v24 * v29 * v36 * v37 * v39 - v12 * v13 * v38 * v40 + v12 * v14 * v38 * v40
          - 3 * v12 * v14 * v15 * v38 * v40 + 2 * v12 * v15 * v16 * v38 * v40
          + v12 * v17 * v38 * v40 - 3 * v12 * v16 * v17 * v38 * v40 - v12 * v18 * v38 * v40
          - v13 * v21 * v38 * v40 + v10 * v13 * v21 * v38 * v40 + v14 * v21 * v38 * v40
          - 3 * v14 * v15 * v21 * v38 * v40 + 3 * v10 * v14 * v15 * v21 * v38 * v40
          + 2 * v15 * v16 * v21 * v38 * v40 - 2 * v10 * v15 * v16 * v21 * v38 * v40
          + v17 * v21 * v38 * v40 - 3 * v16 * v17 * v21 * v38 * v40
          + 3 * v10 * v16 * v17 * v21 * v38 * v40 - v18 * v21 * v38 * v40
          + v10 * v18 * v21 * v38 * v40 + v13 * v29 * v38 * v40 - v14 * v29 * v38 * v40
          + 3 * v14 * v15 * v29 * v38 * v40 - 2 * v15 * v16 * v29 * v38 * v40
          - v17 * v29 * v38 * v40 + 3 * v16 * v17 * v29 * v38 * v40 + v18 * v29 * v38 * v40
          + v14 * v21 * v30 * v38 * v40 + v17 * v21 * v30 * v38 * v40
          - 2 * v10 * v11 * v13 * v38 * v41 + v11 * v12 * v13 * v38 * v41
          + 2 * v10 * v11 * v14 * v38 * v41 - 6 * v10 * v11 * v14 * v15 * v38 * v41
          + 3 * v11 * v12 * v14 * v15 * v38 * v41 + 4 * v10 * v11 * v15 * v16 * v38 * v41
          - 2 * v11 * v12 * v15 * v16 * v38 * v41 + 2 * v10 * v11 * v17 * v38 * v41
          - 6 * v10 * v11 * v16 * v17 * v38 * v41 + 3 * v11 * v12 * v16 * v17 * v38 * v41
          - 2 * v10 * v11 * v18 * v38 * v41 + v11 * v12 * v18 * v38 * v41
          + v10 * v13 * v21 * v38 * v41 + 3 * v10 * v14 * v15 * v21 * v38 * v41
          - 2 * v10 * v15 * v16 * v21 * v38 * v41 + 3 * v10 * v16 * v17 * v21 * v38 * v41
          + v10 * v18 * v21 * v38 * v41 + v14 * v21 * v30 * v38 * v41 + v17 * v21 * v30 * v38 * v41
          + v12 * v14 * v31 * v38 * v41 + v12 * v17 * v31 * v38 * v41
          + 2 * v10 * v11 * v14 * v39 * v42 + v12 * v14 * v39 * v42
          - 3 * v11 * v12 * v14 * v39 * v42 + 2 * v10 * v11 * v17 * v39 * v42
          + v12 * v17 * v39 * v42 - 3 * v11 * v12 * v17 * v39 * v42
          - 2 * v10 * v11 * v16 * v17 * v39 * v42 + 3 * v11 * v12 * v16 * v17 * v39 * v42
          - 2 * v10 * v11 * v18 * v39 * v42 - v12 * v18 * v39 * v42
          + 3 * v11 * v12 * v18 * v39 * v42 + v14 * v21 * v39 * v42
          - 3 * v10 * v14 * v21 * v39 * v42 + v17 * v21 * v39 * v42
          - 3 * v10 * v17 * v21 * v39 * v42 + 3 * v10 * v16 * v17 * v21 * v39 * v42
          - v18 * v21 * v39 * v42 + 3 * v10 * v18 * v21 * v39 * v42 - v14 * v22 * v39 * v42
          - v17 * v22 * v39 * v42 + v16 * v17 * v22 * v39 * v42 + v18 * v22 * v39 * v42
          - v14 * v29 * v39 * v42 - v17 * v29 * v39 * v42 + v16 * v17 * v29 * v39 * v42
          + v18 * v29 * v39 * v42 + v12 * v17 * v32 * v39 * v42 + v17 * v21 * v32 * v39 * v42
          - 2 * v10 * v11 * v14 * v15 * v39 * v43 + 3 * v11 * v12 * v14 * v15 * v39 * v43
          + 4 * v10 * v11 * v15 * v16 * v39 * v43 + 2 * v12 * v15 * v16 * v39 * v43
          - 6 * v11 * v12 * v15 * v16 * v39 * v43 - 2 * v10 * v11 * v16 * v17 * v39 * v43
          + 3 * v11 * v12 * v16 * v17 * v39 * v43 + 3 * v10 * v14 * v15 * v21 * v39 * v43
          + 2 * v15 * v16 * v21 * v39 * v43 - 6 * v10 * v15 * v16 * v21 * v39 * v43
          + 3 * v10 * v16 * v17 * v21 * v39 * v43 + v14 * v15 * v22 * v39 * v43
          - 2 * v15 * v16 * v22 * v39 * v43 + v16 * v17 * v22 * v39 * v43
          + v14 * v15 * v29 * v39 * v43 - 2 * v15 * v16 * v29 * v39 * v43
          + v16 * v17 * v29 * v39 * v43 + v12 * v17 * v32 * v39 * v43 + v17 * v21 * v32 * v39 * v43
          + v12 * v14 * v33 * v39 * v43 + v14 * v21 * v33 * v39 * v43
          - 2 * v10 * v13 * v20 * v34 * v38 * v44 + 2 * v10 * v14 * v20 * v34 * v38 * v44
          - 6 * v10 * v14 * v15 * v20 * v34 * v38 * v44
          + 4 * v10 * v15 * v16 * v20 * v34 * v38 * v44 + 2 * v10 * v17 * v20 * v34 * v38 * v44
          - 6 * v10 * v16 * v17 * v20 * v34 * v38 * v44 - 2 * v10 * v18 * v20 * v34 * v38 * v44
          + 2 * v11 * v13 * v19 * v35 * v38 * v44 - 2 * v11 * v14 * v19 * v35 * v38 * v44
          + 6 * v11 * v14 * v15 * v19 * v35 * v38 * v44
          - 4 * v11 * v15 * v16 * v19 * v35 * v38 * v44 - 2 * v11 * v17 * v19 * v35 * v38 * v44
          + 6 * v11 * v16 * v17 * v19 * v35 * v38 * v44 + 2 * v11 * v18 * v19 * v35 * v38 * v44
          - 2 * v10 * v11 * v13 * v38 * v40 * v44 + 2 * v10 * v11 * v14 * v38 * v40 * v44
          - 6 * v10 * v11 * v14 * v15 * v38 * v40 * v44
          + 4 * v10 * v11 * v15 * v16 * v38 * v40 * v44 + 2 * v10 * v11 * v17 * v38 * v40 * v44
          - 6 * v10 * v11 * v16 * v17 * v38 * v40 * v44 - 2 * v10 * v11 * v18 * v38 * v40 * v44
          + 2 * v10 * v11 * v13 * v38 * v41 * v44 - 2 * v10 * v11 * v14 * v38 * v41 * v44
          + 6 * v10 * v11 * v14 * v15 * v38 * v41 * v44
          - 4 * v10 * v11 * v15 * v16 * v38 * v41 * v44 - 2 * v10 * v11 * v17 * v38 * v41 * v44
          + 6 * v10 * v11 * v16 * v17 * v38 * v41 * v44 + 2 * v10 * v11 * v18 * v38 * v41 * v44
          + 4 * v10 * v11 * v16 * v24 * v36 * v39 * v45 + 2 * v12 * v16 * v24 * v36 * v39 * v45
          - 6 * v11 * v12 * v16 * v24 * v36 * v39 * v45 + 2 * v16 * v21 * v24 * v36 * v39 * v45
          - 6 * v10 * v16 * v21 * v24 * v36 * v39 * v45 - 2 * v16 * v22 * v24 * v36 * v39 * v45
          - 2 * v16 * v24 * v29 * v36 * v39 * v45 - 4 * v10 * v11 * v15 * v23 * v37 * v39 * v45
          - 2 * v12 * v15 * v23 * v37 * v39 * v45 + 6 * v11 * v12 * v15 * v23 * v37 * v39 * v45
          - 2 * v15 * v21 * v23 * v37 * v39 * v45 + 6 * v10 * v15 * v21 * v23 * v37 * v39 * v45
          + 2 * v15 * v22 * v23 * v37 * v39 * v45 + 2 * v15 * v23 * v29 * v37 * v39 * v45
          + 4 * v10 * v11 * v15 * v16 * v39 * v42 * v45 + 2 * v12 * v15 * v16 * v39 * v42 * v45
          - 6 * v11 * v12 * v15 * v16 * v39 * v42 * v45 + 2 * v15 * v16 * v21 * v39 * v42 * v45
          - 6 * v10 * v15 * v16 * v21 * v39 * v42 * v45 - 2 * v15 * v16 * v22 * v39 * v42 * v45
          - 2 * v15 * v16 * v29 * v39 * v42 * v45 - 4 * v10 * v11 * v15 * v16 * v39 * v43 * v45
          - 2 * v12 * v15 * v16 * v39 * v43 * v45 + 6 * v11 * v12 * v15 * v16 * v39 * v43 * v45
          - 2 * v15 * v16 * v21 * v39 * v43 * v45 + 6 * v10 * v15 * v16 * v21 * v39 * v43 * v45
          + 2 * v15 * v16 * v22 * v39 * v43 * v45 + 2 * v15 * v16 * v29 * v39 * v43 * v45
          + 2 * v13 * v20 * v34 * v35 * v38 * v6 - 2 * v14 * v20 * v34 * v35 * v38 * v6
          + 6 * v14 * v15 * v20 * v34 * v35 * v38 * v6 - 4 * v15 * v16 * v20 * v34 * v35 * v38 * v6
          - 2 * v17 * v20 * v34 * v35 * v38 * v6 + 6 * v16 * v17 * v20 * v34 * v35 * v38 * v6
          + 2 * v18 * v20 * v34 * v35 * v38 * v6 - 2 * v13 * v28 * v34 * v35 * v38 * v6
          + 2 * v14 * v28 * v34 * v35 * v38 * v6 - 6 * v14 * v15 * v28 * v34 * v35 * v38 * v6
          + 4 * v15 * v16 * v28 * v34 * v35 * v38 * v6 + 2 * v17 * v28 * v34 * v35 * v38 * v6
          - 6 * v16 * v17 * v28 * v34 * v35 * v38 * v6 - 2 * v18 * v28 * v34 * v35 * v38 * v6
          + 2 * v13 * v21 * v35 * v38 * v44 * v6 - 2 * v14 * v21 * v35 * v38 * v44 * v6
          + 6 * v14 * v15 * v21 * v35 * v38 * v44 * v6 - 4 * v15 * v16 * v21 * v35 * v38 * v44 * v6
          - 2 * v17 * v21 * v35 * v38 * v44 * v6 + 6 * v16 * v17 * v21 * v35 * v38 * v44 * v6
          + 2 * v18 * v21 * v35 * v38 * v44 * v6 - 2 * v13 * v20 * v34 * v35 * v38 * v44 * v6
          + 2 * v14 * v20 * v34 * v35 * v38 * v44 * v6
          - 6 * v14 * v15 * v20 * v34 * v35 * v38 * v44 * v6
          + 4 * v15 * v16 * v20 * v34 * v35 * v38 * v44 * v6
          + 2 * v17 * v20 * v34 * v35 * v38 * v44 * v6
          - 6 * v16 * v17 * v20 * v34 * v35 * v38 * v44 * v6
          - 2 * v18 * v20 * v34 * v35 * v38 * v44 * v6 + 2 * v12 * v13 * v34 * v38 * v7
          - 2 * v12 * v14 * v34 * v38 * v7 + 6 * v12 * v14 * v15 * v34 * v38 * v7
          - 4 * v12 * v15 * v16 * v34 * v38 * v7 - 2 * v12 * v17 * v34 * v38 * v7
          + 6 * v12 * v16 * v17 * v34 * v38 * v7 + 2 * v12 * v18 * v34 * v38 * v7
          - 2 * v13 * v19 * v34 * v35 * v38 * v7 + 2 * v14 * v19 * v34 * v35 * v38 * v7
          - 6 * v14 * v15 * v19 * v34 * v35 * v38 * v7 + 4 * v15 * v16 * v19 * v34 * v35 * v38 * v7
          + 2 * v17 * v19 * v34 * v35 * v38 * v7 - 6 * v16 * v17 * v19 * v34 * v35 * v38 * v7
          - 2 * v18 * v19 * v34 * v35 * v38 * v7 - 2 * v12 * v13 * v34 * v38 * v44 * v7
          + 2 * v12 * v14 * v34 * v38 * v44 * v7 - 6 * v12 * v14 * v15 * v34 * v38 * v44 * v7
          + 4 * v12 * v15 * v16 * v34 * v38 * v44 * v7 + 2 * v12 * v17 * v34 * v38 * v44 * v7
          - 6 * v12 * v16 * v17 * v34 * v38 * v44 * v7 - 2 * v12 * v18 * v34 * v38 * v44 * v7
          + 2 * v13 * v19 * v34 * v35 * v38 * v44 * v7 - 2 * v14 * v19 * v34 * v35 * v38 * v44 * v7
          + 6 * v14 * v15 * v19 * v34 * v35 * v38 * v44 * v7
          - 4 * v15 * v16 * v19 * v34 * v35 * v38 * v44 * v7
          - 2 * v17 * v19 * v34 * v35 * v38 * v44 * v7
          + 6 * v16 * v17 * v19 * v34 * v35 * v38 * v44 * v7
          + 2 * v18 * v19 * v34 * v35 * v38 * v44 * v7 - 4 * v10 * v11 * v24 * v36 * v37 * v39 * v8
          - 2 * v12 * v24 * v36 * v37 * v39 * v8 + 6 * v11 * v12 * v24 * v36 * v37 * v39 * v8
          - 2 * v21 * v24 * v36 * v37 * v39 * v8 + 6 * v10 * v21 * v24 * v36 * v37 * v39 * v8
          + 2 * v22 * v24 * v36 * v37 * v39 * v8 + 4 * v10 * v11 * v26 * v36 * v37 * v39 * v8
          + 2 * v12 * v26 * v36 * v37 * v39 * v8 - 6 * v11 * v12 * v26 * v36 * v37 * v39 * v8
          + 2 * v21 * v26 * v36 * v37 * v39 * v8 - 6 * v10 * v21 * v26 * v36 * v37 * v39 * v8
          - 2 * v22 * v26 * v36 * v37 * v39 * v8 + 2 * v24 * v29 * v36 * v37 * v39 * v8
          - 2 * v26 * v29 * v36 * v37 * v39 * v8 - 4 * v10 * v11 * v17 * v37 * v39 * v45 * v8
          - 2 * v12 * v17 * v37 * v39 * v45 * v8 + 6 * v11 * v12 * v17 * v37 * v39 * v45 * v8
          - 2 * v17 * v21 * v37 * v39 * v45 * v8 + 6 * v10 * v17 * v21 * v37 * v39 * v45 * v8
          + 2 * v17 * v22 * v37 * v39 * v45 * v8 + 2 * v17 * v29 * v37 * v39 * v45 * v8
          + 4 * v10 * v11 * v24 * v36 * v37 * v39 * v45 * v8
          + 2 * v12 * v24 * v36 * v37 * v39 * v45 * v8
          - 6 * v11 * v12 * v24 * v36 * v37 * v39 * v45 * v8
          + 2 * v21 * v24 * v36 * v37 * v39 * v45 * v8
          - 6 * v10 * v21 * v24 * v36 * v37 * v39 * v45 * v8
          - 2 * v22 * v24 * v36 * v37 * v39 * v45 * v8 - 2 * v24 * v29 * v36 * v37 * v39 * v45 * v8
          - 4 * v10 * v11 * v14 * v36 * v39 * v9 - 2 * v12 * v14 * v36 * v39 * v9
          + 6 * v11 * v12 * v14 * v36 * v39 * v9 - 2 * v14 * v21 * v36 * v39 * v9
          + 6 * v10 * v14 * v21 * v36 * v39 * v9 + 2 * v14 * v22 * v36 * v39 * v9
          + 2 * v14 * v29 * v36 * v39 * v9 + 4 * v10 * v11 * v23 * v36 * v37 * v39 * v9
          + 2 * v12 * v23 * v36 * v37 * v39 * v9 - 6 * v11 * v12 * v23 * v36 * v37 * v39 * v9
          + 2 * v21 * v23 * v36 * v37 * v39 * v9 - 6 * v10 * v21 * v23 * v36 * v37 * v39 * v9
          - 2 * v22 * v23 * v36 * v37 * v39 * v9 - 2 * v23 * v29 * v36 * v37 * v39 * v9
          + 4 * v10 * v11 * v14 * v36 * v39 * v45 * v9 + 2 * v12 * v14 * v36 * v39 * v45 * v9
          - 6 * v11 * v12 * v14 * v36 * v39 * v45 * v9 + 2 * v14 * v21 * v36 * v39 * v45 * v9
          - 6 * v10 * v14 * v21 * v36 * v39 * v45 * v9 - 2 * v14 * v22 * v36 * v39 * v45 * v9
          - 2 * v14 * v29 * v36 * v39 * v45 * v9 - 4 * v10 * v11 * v23 * v36 * v37 * v39 * v45 * v9
          - 2 * v12 * v23 * v36 * v37 * v39 * v45 * v9
          + 6 * v11 * v12 * v23 * v36 * v37 * v39 * v45 * v9
          - 2 * v21 * v23 * v36 * v37 * v39 * v45 * v9
          + 6 * v10 * v21 * v23 * v36 * v37 * v39 * v45 * v9
          + 2 * v22 * v23 * v36 * v37 * v39 * v45 * v9 + 2 * v23 * v29 * v36 * v37 * v39 * v45 * v9)
         / ((-1 + v10 + v11) * pow(v10 + v11, 2) * (-1 + v15 + v16) * pow(v15 + v16, 2));
}