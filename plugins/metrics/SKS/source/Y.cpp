#include "SKS.hpp"

auto grlensing::SKS::Y(unsigned bhIdx, double t, double x, double y, double z) const noexcept
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

  return (v1 * v2 * v3 + v0 * v7 + v0 * v6 * sqrt(1 - v4 + v7) - v1 * v2 * v3 * v8 - v2 * v3 * x
          + v2 * v3 * v8 * x + v5 * y + v4 * v8 * y)
         / (v4 + v5);
}