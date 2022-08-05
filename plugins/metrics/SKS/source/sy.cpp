#include "SKS.hpp"

auto grlensing::SKS::sy(unsigned bhIdx, double t) const noexcept -> double {

  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  return -(pow(-1, bhIdx) * b * sin(sqrt((M1 + M2) / pow(b, 3)) * t)) / 2.;
}