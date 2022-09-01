#include "SKS.hpp"

auto grlensing::SKS::dX_dt(unsigned bhIdx, double t, double x, double y, double z) const noexcept
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
  const double v8{pow(v4, 3)};
  const double v9{pow(v3, 3)};
  const double v10{pow(v4, 4)};
  const double v11{pow(v4, 5)};
  const double v12{pow(v3, 5)};
  const double v13{pow(v3, 4)};
  const double v14{-v6};
  const double v15{-v7};
  const double v16{1 + v14 + v15};
  const double v17{sqrt(v16)};

  return (-v11 + v0 * v11 * v2 - v1 * v10 * v2 * v3 - v13 * v4 + pow(v3, 6) * v4 + pow(v4, 7)
          + v1 * v11 * v5 + v0 * v12 * v5 + 2 * v1 * v13 * v4 * v5
          + v0 * v14 * v2 * v4 * sqrt(1 + v15 - v6) + 3 * v11 * v6 + v0 * v2 * v4 * v6
          - 2 * v1 * v4 * v5 * v6 + 2 * v1 * v17 * v4 * v5 * v6
          + v0 * v15 * v3 * v5 * sqrt(1 + v14 - v7) + 2 * v1 * v2 * v3 * v7
          - 2 * v1 * v17 * v2 * v3 * v7 + v0 * v3 * v5 * v7 + 3 * v13 * v8 - v0 * v2 * v8
          + v0 * v17 * v2 * v8 - 2 * v6 * v8 + v0 * v2 * v6 * v8 + 3 * v1 * v5 * v6 * v8
          + v1 * v15 * v2 * v9 - v0 * v5 * v9 + v0 * v17 * v5 * v9 + v0 * v5 * v7 * v9
          + v10 * v2 * v3 * x - v11 * v5 * x - 2 * v13 * v4 * v5 * x + 2 * v4 * v5 * v6 * x
          - 2 * v17 * v4 * v5 * v6 * x - 2 * v2 * v3 * v7 * x + 2 * v17 * v2 * v3 * v7 * x
          - 3 * v5 * v6 * v8 * x + v2 * v7 * v9 * x - v11 * v2 * y + v14 * v2 * v4 * y
          - v12 * v5 * y + v15 * v3 * v5 * y + v17 * v2 * v4 * v6 * y + v17 * v3 * v5 * v7 * y
          + v2 * v8 * y + v14 * v2 * v8 * y - v17 * v2 * v8 * y + v5 * v9 * y + v15 * v5 * v9 * y
          - v17 * v5 * v9 * y)
         / (sqrt(v16) * pow(v6 + v7, 2));
}