#include "SKS.hpp"

auto grlensing::SKS::T(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0 = dsy_dt(bhIdx, t);
  const double v1 = dsx_dt(bhIdx, t);

  return sqrt(1 - pow(v0, 2) - pow(v1, 2))
         * (t + v1 * x + v0 * y - v1 * sx(bhIdx, t) - v0 * sy(bhIdx, t));
}