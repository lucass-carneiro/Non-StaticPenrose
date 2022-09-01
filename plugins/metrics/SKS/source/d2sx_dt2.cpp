#include "SKS.hpp"

auto grlensing::SKS::d2sx_dt2(unsigned bhIdx, double t) const noexcept -> double {

  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  const double v0{M1 + M2};

  return (pow(-1, bhIdx) * v0 * cos(t * sqrt(v0 / pow(b, 3)))) / (2. * pow(b, 2));
}