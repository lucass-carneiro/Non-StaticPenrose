#include "SKS.hpp"

auto grlensing::SKS::dsy_dt(unsigned bhIdx, double t) const noexcept -> double {

  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  const double v0{sqrt((M1 + M2) / pow(b, 3))};

  return -(pow(-1, bhIdx) * b * v0 * cos(t * v0)) / 2.;
}