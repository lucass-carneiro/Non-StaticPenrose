#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::lapse(double t, double x, double y, double z) -> double {
  using namespace sks_aux;
  using std::sqrt;

  const auto llg = llg_SKS(M1, M2, a1, a2, b, t, x, y, z);
  const auto uug = inverse_symmetric_4x4(llg);

  return sqrt(-1.0 / uug[0][0]);
}