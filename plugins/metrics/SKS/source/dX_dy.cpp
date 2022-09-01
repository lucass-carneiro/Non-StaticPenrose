#include "SKS.hpp"

auto grlensing::SKS::dX_dy(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0{dsy_dt(bhIdx, t)};
  const double v1{dsx_dt(bhIdx, t)};
  const double v2{pow(v0, 2)};
  const double v3{pow(v1, 2)};

  return (v0 * v1 * (-1 + sqrt(1 - v2 - v3))) / (v2 + v3);
}