#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::dllgSKS_03_dx(double t, double x, double y, double z) const noexcept
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
  const double v10{dX_dx(2, t, x, y, z)};
  const double v11{dX_dx(1, t, x, y, z)};
  const double v12{dZ_dx(2, t, x, y, z)};
  const double v13{dZ_dx(1, t, x, y, z)};
  const double v14{dY_dx(2, t, x, y, z)};
  const double v15{dY_dx(1, t, x, y, z)};
  const double v16{l3_KS(a2, v2, v1, v0)};
  const double v17{l2_KS(a2, v2, v1, v0)};
  const double v18{l1_KS(a2, v2, v1, v0)};
  const double v19{l3_KS(a1, v5, v4, v3)};
  const double v20{l2_KS(a1, v5, v4, v3)};
  const double v21{l1_KS(a1, v5, v4, v3)};
  const double v22{H_KS(M2, a2, v2, v1, v0)};
  const double v23{H_KS(M1, a1, v5, v4, v3)};
  const double v24{dl3_KS_dZ(a2, v2, v1, v0)};
  const double v25{dl3_KS_dZ(a1, v5, v4, v3)};
  const double v26{dl3_KS_dY(a2, v2, v1, v0)};
  const double v27{dl3_KS_dY(a1, v5, v4, v3)};
  const double v28{dl3_KS_dX(a2, v2, v1, v0)};
  const double v29{dl3_KS_dX(a1, v5, v4, v3)};
  const double v30{1 - pow(v7, 2) - pow(v8, 2)};
  const double v31{1 - pow(v6, 2) - pow(v9, 2)};
  const double v32{dH_KS_dX(M2, a2, v2, v1, v0)};
  const double v33{dH_KS_dX(M1, a1, v5, v4, v3)};
  const double v34{dH_KS_dZ(M2, a2, v2, v1, v0)};
  const double v35{dH_KS_dZ(M1, a1, v5, v4, v3)};
  const double v36{dH_KS_dY(M2, a2, v2, v1, v0)};
  const double v37{dH_KS_dY(M1, a1, v5, v4, v3)};
  const double v38{sqrt(v30)};
  const double v39{sqrt(v31)};

  return (v12 * v22 * v24 * v38 + v14 * v22 * v26 * v38 + v10 * v22 * v28 * v38
          + v10 * v16 * v32 * v38 + v12 * v16 * v34 * v38 + v14 * v16 * v36 * v38
          + v13 * v23 * v25 * v39 + v15 * v23 * v27 * v39 + v11 * v23 * v29 * v39
          + v11 * v19 * v33 * v39 + v13 * v19 * v35 * v39 + v15 * v19 * v37 * v39
          - v12 * v17 * v22 * v24 * v38 * v6 - v14 * v17 * v22 * v26 * v38 * v6
          - v10 * v17 * v22 * v28 * v38 * v6 - v10 * v16 * v17 * v32 * v38 * v6
          - v12 * v16 * v17 * v34 * v38 * v6 - v14 * v16 * v17 * v36 * v38 * v6
          - v13 * v20 * v23 * v25 * v39 * v7 - v15 * v20 * v23 * v27 * v39 * v7
          - v11 * v20 * v23 * v29 * v39 * v7 - v11 * v19 * v20 * v33 * v39 * v7
          - v13 * v19 * v20 * v35 * v39 * v7 - v15 * v19 * v20 * v37 * v39 * v7
          - v13 * v21 * v23 * v25 * v39 * v8 - v15 * v21 * v23 * v27 * v39 * v8
          - v11 * v21 * v23 * v29 * v39 * v8 - v11 * v19 * v21 * v33 * v39 * v8
          - v13 * v19 * v21 * v35 * v39 * v8 - v15 * v19 * v21 * v37 * v39 * v8
          - v12 * v18 * v22 * v24 * v38 * v9 - v14 * v18 * v22 * v26 * v38 * v9
          - v10 * v18 * v22 * v28 * v38 * v9 - v10 * v16 * v18 * v32 * v38 * v9
          - v12 * v16 * v18 * v34 * v38 * v9 - v14 * v16 * v18 * v36 * v38 * v9
          - v13 * v19 * v23 * v39 * v8 * dl1_KS_dZ(a1, v5, v4, v3)
          - v12 * v16 * v22 * v38 * v9 * dl1_KS_dZ(a2, v2, v1, v0)
          - v13 * v19 * v23 * v39 * v7 * dl2_KS_dZ(a1, v5, v4, v3)
          - v12 * v16 * v22 * v38 * v6 * dl2_KS_dZ(a2, v2, v1, v0)
          - v15 * v19 * v23 * v39 * v8 * dl1_KS_dY(a1, v5, v4, v3)
          - v14 * v16 * v22 * v38 * v9 * dl1_KS_dY(a2, v2, v1, v0)
          - v15 * v19 * v23 * v39 * v7 * dl2_KS_dY(a1, v5, v4, v3)
          - v14 * v16 * v22 * v38 * v6 * dl2_KS_dY(a2, v2, v1, v0)
          - v11 * v19 * v23 * v39 * v8 * dl1_KS_dX(a1, v5, v4, v3)
          - v10 * v16 * v22 * v38 * v9 * dl1_KS_dX(a2, v2, v1, v0)
          - v11 * v19 * v23 * v39 * v7 * dl2_KS_dX(a1, v5, v4, v3)
          - v10 * v16 * v22 * v38 * v6 * dl2_KS_dX(a2, v2, v1, v0))
         / (sqrt(v30) * sqrt(v31));
}