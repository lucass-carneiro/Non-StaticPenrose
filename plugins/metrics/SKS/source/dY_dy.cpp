#include "SKS.hpp"

auto grlensing::SKS::dY_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0 = pow(dsy_dt(bhIdx, t), 2);
  const double v1 = pow(dsx_dt(bhIdx, t), 2);

  return (v0 * sqrt(1 - v0 - v1) + v1) / (v0 + v1);
}