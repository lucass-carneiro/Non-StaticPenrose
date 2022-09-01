#include "SKS.hpp"

auto grlensing::SKS::dT_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0{dsy_dt(bhIdx, t)};

  return v0 * sqrt(1 - pow(v0, 2) - pow(dsx_dt(bhIdx, t), 2));
}