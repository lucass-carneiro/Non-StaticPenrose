#include "SKS.hpp"
#include "aux_functions.hpp"

auto grlensing::SKS::llgSKS_33(double t, double x, double y, double z) const noexcept -> double {

  using namespace sks_aux;
  using std::pow;
  using std::sqrt;

  const double v0 = Z(2, t, x, y, z);
  const double v1 = Y(2, t, x, y, z);
  const double v2 = X(2, t, x, y, z);
  const double v3 = Z(1, t, x, y, z);
  const double v4 = Y(1, t, x, y, z);
  const double v5 = X(1, t, x, y, z);

  return 1 + H_KS(M1, a1, v5, v4, v3) * pow(l3_KS(a1, v5, v4, v3), 2)
         + H_KS(M2, a2, v2, v1, v0) * pow(l3_KS(a2, v2, v1, v0), 2);
}