#include "SKS.hpp"

auto grlensing::SKS::dT_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0{sy(bhIdx, t)};
  const double v1{sx(bhIdx, t)};
  const double v2{d2sy_dt2(bhIdx, t)};
  const double v3{dsy_dt(bhIdx, t)};
  const double v4{dsx_dt(bhIdx, t)};
  const double v5{d2sx_dt2(bhIdx, t)};
  const double v6{pow(v3, 2)};
  const double v7{pow(v4, 2)};

  return (1 - v0 * v2 - t * v2 * v3 + pow(v3, 4) + v1 * v2 * v3 * v4 + pow(v4, 4) - v1 * v5
          - t * v4 * v5 + v0 * v3 * v4 * v5 - 2 * v6 + 2 * v0 * v2 * v6 + v1 * v5 * v6 - 2 * v7
          + v0 * v2 * v7 + 2 * v1 * v5 * v7 + 2 * v6 * v7 - v2 * v3 * v4 * x + v5 * x - v5 * v6 * x
          - 2 * v5 * v7 * x + v2 * y - v3 * v4 * v5 * y - 2 * v2 * v6 * y - v2 * v7 * y)
         / sqrt(1 - v6 - v7);
}