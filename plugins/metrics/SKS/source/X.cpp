#include "SKS.hpp"

auto grlensing::SKS::X(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0{sy(bhIdx, t)};
  const double v1{sx(bhIdx, t)};
  const double v2{dsy_dt(bhIdx, t)};
  const double v3{dsx_dt(bhIdx, t)};
  const double v4{pow(v2, 2)};
  const double v5{pow(v3, 2)};
  const double v6{-v4};
  const double v7{-v5};
  const double v8{sqrt(1 + v6 + v7)};

  return (v0 * v2 * v3 + v1 * v6 + v1 * sqrt(1 - v5 + v6) * v7 - v0 * v2 * v3 * v8 + v4 * x
          + v5 * v8 * x - v2 * v3 * y + v2 * v3 * v8 * y)
         / (v4 + v5);
}