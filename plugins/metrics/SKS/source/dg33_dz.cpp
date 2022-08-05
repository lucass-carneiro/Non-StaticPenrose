#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::dllgSKS_33_dz(double t, double x, double y, double z) const noexcept
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
  const double v6 = dX_dz(2, t, x, y, z);
  const double v7 = dX_dz(1, t, x, y, z);
  const double v8 = dY_dz(2, t, x, y, z);
  const double v9 = dY_dz(1, t, x, y, z);
  const double v10 = dZ_dz(2, t, x, y, z);
  const double v11 = dZ_dz(1, t, x, y, z);
  const double v12 = l3_KS(a2, v2, v1, v0);
  const double v13 = l3_KS(a1, v5, v4, v3);
  const double v14 = H_KS(M2, a2, v2, v1, v0);
  const double v15 = H_KS(M1, a1, v5, v4, v3);
  const double v16 = pow(v12, 2);
  const double v17 = pow(v13, 2);

  return 2 * v11 * v13 * v15 * dl2_KS_dZ(a1, v5, v4, v3)
         + 2 * v10 * v12 * v14 * dl2_KS_dZ(a2, v2, v1, v0)
         + 2 * v13 * v15 * v9 * dl2_KS_dY(a1, v5, v4, v3)
         + 2 * v12 * v14 * v8 * dl2_KS_dY(a2, v2, v1, v0)
         + 2 * v13 * v15 * v7 * dl2_KS_dX(a1, v5, v4, v3)
         + 2 * v12 * v14 * v6 * dl2_KS_dX(a2, v2, v1, v0) + v11 * v17 * dH_KS_dZ(M1, a1, v5, v4, v3)
         + v10 * v16 * dH_KS_dZ(M2, a2, v2, v1, v0) + v17 * v9 * dH_KS_dY(M1, a1, v5, v4, v3)
         + v16 * v8 * dH_KS_dY(M2, a2, v2, v1, v0) + v17 * v7 * dH_KS_dX(M1, a1, v5, v4, v3)
         + v16 * v6 * dH_KS_dX(M2, a2, v2, v1, v0);
}