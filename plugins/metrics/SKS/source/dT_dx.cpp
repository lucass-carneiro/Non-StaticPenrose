#include "SKS.hpp"

auto grlensing::SKS::dT_dx(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0 = dsx_dt(bhIdx, t);

  return v0 * sqrt(1 - pow(v0, 2) - pow(dsy_dt(bhIdx, t), 2));
}