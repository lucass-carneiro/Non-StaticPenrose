#include "SKS.hpp"

auto grlensing::SKS::dY_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
    -> double {

  using std::pow;
  using std::sqrt;

  const double v0{sx(bhIdx, t)};
  const double v1{sy(bhIdx, t)};
  const double v2{d2sy_dt2(bhIdx, t)};
  const double v3{dsy_dt(bhIdx, t)};
  const double v4{dsx_dt(bhIdx, t)};
  const double v5{d2sx_dt2(bhIdx, t)};
  const double v6{pow(v3, 2)};
  const double v7{pow(v4, 2)};
  const double v8{pow(v4, 3)};
  const double v9{pow(v3, 5)};
  const double v10{pow(v3, 3)};
  const double v11{pow(v4, 4)};
  const double v12{pow(v4, 5)};
  const double v13{pow(v3, 4)};
  const double v14{-v6};
  const double v15{-v7};
  const double v16{1 + v14 + v15};
  const double v17{sqrt(v16)};

  return (3 * v10 * v11 + v0 * v12 * v2 - v11 * v3 + 2 * v1 * v11 * v2 * v3 + pow(v3, 7)
          + v3 * pow(v4, 6) - v0 * v10 * v5 + v0 * v10 * v17 * v5 - v1 * v13 * v4 * v5
          + v0 * v14 * v2 * v4 * sqrt(1 + v15 - v6) + v0 * v2 * v4 * v6 + 2 * v1 * v4 * v5 * v6
          - 2 * v1 * v17 * v4 * v5 * v6 + v0 * v15 * v3 * v5 * sqrt(1 + v14 - v7) - 2 * v10 * v7
          + 3 * v1 * v10 * v2 * v7 - 2 * v1 * v2 * v3 * v7 + 2 * v1 * v17 * v2 * v3 * v7
          + v0 * v10 * v5 * v7 + v0 * v3 * v5 * v7 - v0 * v2 * v8 + v0 * v17 * v2 * v8
          + v1 * v14 * v5 * v8 + v0 * v2 * v6 * v8 - v9 + v1 * v2 * v9 + v0 * v5 * v9 + 3 * v7 * v9
          - v12 * v2 * x + v14 * v2 * v4 * x + v10 * v5 * x + v10 * v15 * v5 * x
          - v10 * v17 * v5 * x + v15 * v3 * v5 * x + v17 * v2 * v4 * v6 * x + v17 * v3 * v5 * v7 * x
          + v2 * v8 * x + v14 * v2 * v8 * x - v17 * v2 * v8 * x - v5 * v9 * x
          - 2 * v11 * v2 * v3 * y + v13 * v4 * v5 * y - 2 * v4 * v5 * v6 * y
          + 2 * v17 * v4 * v5 * v6 * y - 3 * v10 * v2 * v7 * y + 2 * v2 * v3 * v7 * y
          - 2 * v17 * v2 * v3 * v7 * y + v5 * v6 * v8 * y - v2 * v9 * y)
         / (sqrt(v16) * pow(v6 + v7, 2));
}