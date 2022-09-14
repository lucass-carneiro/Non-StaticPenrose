#ifndef GRLENSING_SKS_AUX_FUNCTIONS
#define GRLENSING_SKS_AUX_FUNCTIONS

#include "../../../../include/tensor_types.hpp"

#include <cmath>
#include <cstddef>

namespace sks_aux {

inline auto H_KS(double M, double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};
  const double v5{v0 * v3};
  const double v6{v0 + v1 + v2 + v4};
  const double v7{sqrt(v5 + pow(v6, 2) / 4.)};

  return (sqrt(2) * M * pow(sqrt(pow(v0 + v1 + v2 + v4, 2) + 4 * v5) + v6, 1.5))
         / (pow(a, 4) + 2 * v1 * v2 - 2 * v1 * v3 - 2 * v2 * v3 + 2 * v5 + 2 * (v1 + v2 - v3) * v7
            + 2 * v0 * (v1 + v2 + v7) + pow(X, 4) + pow(Y, 4) + pow(Z, 4));
}

inline auto l1_KS(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{v0 + v1 + v2 - v3};
  const double v5{sqrt(v0 * v3 + pow(v4, 2) / 4.)};

  return (2 * sqrt(v4 / 2. + v5) * X + 2 * a * Y) / (v0 + v1 + v2 + v3 + 2 * v5);
}

inline auto l2_KS(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{v0 + v1 + v2 - v3};
  const double v5{sqrt(v0 * v3 + pow(v4, 2) / 4.)};

  return (-2 * a * X + 2 * sqrt(v4 / 2. + v5) * Y) / (v0 + v1 + v2 + v3 + 2 * v5);
}

inline auto l3_KS(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};

  return (sqrt(2) * Z)
         / sqrt(v0 + v1 + v2 + v4 + 2 * sqrt(v0 * v3 + pow(v0 + v1 + v2 + v4, 2) / 4.));
}

inline auto dH_KS_dX(double M, double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Z, 4)};
  const double v5{pow(Y, 4)};
  const double v6{pow(X, 4)};
  const double v7{pow(a, 4)};
  const double v8{-v3};
  const double v9{v0 * v3};
  const double v10{2 * v0 * v1};
  const double v11{2 * v0 * v2};
  const double v12{2 * v9};
  const double v13{2 * v1 * v2};
  const double v14{-2 * v1 * v3};
  const double v15{-2 * v2 * v3};
  const double v16{v0 + v1 + v2 + v8};
  const double v17{sqrt(pow(v16, 2) / 4. + v9)};
  const double v18{2 * v0 * v17};
  const double v19{2 * v1 * v17};
  const double v20{2 * v17 * v2};
  const double v21{-2 * v17 * v3};

  return -(
      (sqrt(2) * M
       * (v10 + v11 + v13 + v14 + v15 + v18 + v19 + v20 + v21 + v4 + v5 + v6 + v7 - 6 * v9)
       * pow(v16 + sqrt(pow(v0 + v1 + v2 + v8, 2) + 4 * v9), 1.5) * X)
      / (sqrt(v10 + v11 + v12 + v13 + v14 + v15 + v4 + v5 + v6 + v7)
         * pow(v10 + v11 + v12 + v13 + v14 + v15 + v18 + v19 + v20 + v21 + v4 + v5 + v6 + v7, 2)));
}

inline auto dl1_KS_dX(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Z, 4)};
  const double v5{pow(Y, 4)};
  const double v6{pow(a, 4)};
  const double v7{-v3};
  const double v8{sqrt(2)};
  const double v9{v0 * v3};
  const double v10{v0 + v1 + v2 + v7};
  const double v11{sqrt(pow(v10, 2) / 4. + v9)};
  const double v12{2 * v11};
  const double v13{sqrt(v10 + sqrt(pow(v0 + v1 + v2 + v7, 2) + 4 * v9))};

  return (v13
              * (v0 * (2 * v1 + v12 + v2) + v1 * (v12 + v2 - 2 * v3) + (v12 + v2) * v3 + v4 + v5
                 + v6)
              * v8
          + 2 * sqrt(v10 + 2 * sqrt(v0 * v3 + pow(v0 + v1 + v2 + v7, 2) / 4.)) * v8 * v9
          - 4 * a * X * Y * (-pow(a, 2) + v0 + 2 * v11 + pow(X, 2) + pow(Y, 2)))
         / (pow(v0 + v1 + v12 + v2 + v3, 2)
            * sqrt(2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v4 + v5 + v6 + 2 * v9
                   + pow(X, 4)));
}

inline auto dl2_KS_dX(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(a, 3)};
  const double v5{pow(Z, 4)};
  const double v6{pow(Y, 4)};
  const double v7{pow(X, 4)};
  const double v8{-v3};
  const double v9{sqrt(2)};
  const double v10{v0 * v3};
  const double v11{v0 + v1 + v2 + v8};
  const double v12{sqrt(v10 + pow(v11, 2) / 4.)};
  const double v13{sqrt(v11 + sqrt(4 * v10 + pow(v0 + v1 + v2 + v8, 2)))};

  return (-2
              * (pow(a, 5) + 2 * (v0 - v1 + v12) * v4
                 + a * (2 * v1 * v12 + 2 * v0 * (v1 + v12) - 2 * v12 * v2 + v5 + v6 - v7))
          - v13 * v9 * X * (v0 + 2 * v12 - 3 * v3 + pow(X, 2)) * Y - v13 * v9 * X * pow(Y, 3))
         / (pow(v0 + v1 + 2 * v12 + v2 + v3, 2)
            * sqrt(pow(a, 4) + 2 * v10 + 2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v5
                   + v6 + v7));
}

inline auto dl3_KS_dX(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};

  return -((sqrt(2) * X * Z)
           / (sqrt(v0 + v1 + v2 + v4 + 2 * sqrt(v0 * v3 + pow(v0 + v1 + v2 + v4, 2) / 4.))
              * sqrt(pow(a, 4) + 2 * v1 * v2 - 2 * (v1 + v2) * v3 + 2 * v0 * (v1 + v2 + v3)
                     + pow(X, 4) + pow(Y, 4) + pow(Z, 4))));
}

inline auto dH_KS_dY(double M, double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Z, 4)};
  const double v5{pow(Y, 4)};
  const double v6{pow(X, 4)};
  const double v7{pow(a, 4)};
  const double v8{-v3};
  const double v9{v0 * v3};
  const double v10{2 * v0 * v1};
  const double v11{2 * v0 * v2};
  const double v12{2 * v9};
  const double v13{2 * v1 * v2};
  const double v14{-2 * v1 * v3};
  const double v15{-2 * v2 * v3};
  const double v16{v0 + v1 + v2 + v8};
  const double v17{sqrt(pow(v16, 2) / 4. + v9)};
  const double v18{2 * v0 * v17};
  const double v19{2 * v1 * v17};
  const double v20{2 * v17 * v2};
  const double v21{-2 * v17 * v3};

  return -(
      (sqrt(2) * M
       * (v10 + v11 + v13 + v14 + v15 + v18 + v19 + v20 + v21 + v4 + v5 + v6 + v7 - 6 * v9)
       * pow(v16 + sqrt(pow(v0 + v1 + v2 + v8, 2) + 4 * v9), 1.5) * Y)
      / (sqrt(v10 + v11 + v12 + v13 + v14 + v15 + v4 + v5 + v6 + v7)
         * pow(v10 + v11 + v12 + v13 + v14 + v15 + v18 + v19 + v20 + v21 + v4 + v5 + v6 + v7, 2)));
}

inline auto dl1_KS_dY(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(a, 3)};
  const double v5{pow(Z, 4)};
  const double v6{pow(Y, 4)};
  const double v7{pow(X, 4)};
  const double v8{-v3};
  const double v9{sqrt(2)};
  const double v10{v0 * v3};
  const double v11{v0 + v1 + v2 + v8};
  const double v12{sqrt(v10 + pow(v11, 2) / 4.)};
  const double v13{sqrt(v11 + sqrt(4 * v10 + pow(v0 + v1 + v2 + v8, 2)))};

  return (2 * pow(a, 5) + 4 * (v0 + v12 - v2) * v4
          + 2 * a * (-2 * v1 * v12 + 2 * v12 * v2 + 2 * v0 * (v12 + v2) + v5 - v6 + v7)
          - v13 * v9 * X * (v0 + 2 * v12 - 3 * v3 + pow(X, 2)) * Y - v13 * v9 * X * pow(Y, 3))
         / (pow(v0 + v1 + 2 * v12 + v2 + v3, 2)
            * sqrt(pow(a, 4) + 2 * v10 + 2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v5
                   + v6 + v7));
}

inline auto dl2_KS_dY(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Z, 4)};
  const double v5{pow(X, 4)};
  const double v6{pow(a, 4)};
  const double v7{-v3};
  const double v8{sqrt(2)};
  const double v9{v0 * v3};
  const double v10{v0 + v1 + v2 + v7};
  const double v11{sqrt(pow(v10, 2) / 4. + v9)};
  const double v12{2 * v11};
  const double v13{sqrt(v10 + sqrt(pow(v0 + v1 + v2 + v7, 2) + 4 * v9))};

  return (v13
              * ((v1 + v12) * v2 + v0 * (v1 + v12 + 2 * v2) + (v1 + v12 - 2 * v2) * v3 + v4 + v5
                 + v6)
              * v8
          + 2 * sqrt(v10 + 2 * sqrt(v0 * v3 + pow(v0 + v1 + v2 + v7, 2) / 4.)) * v8 * v9
          + 4 * a * X * Y * (-pow(a, 2) + v0 + 2 * v11 + pow(X, 2) + pow(Y, 2)))
         / (pow(v0 + v1 + v12 + v2 + v3, 2)
            * sqrt(2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v4 + v5 + v6 + 2 * v9
                   + pow(Y, 4)));
}

inline auto dl3_KS_dY(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};

  return -((sqrt(2) * Y * Z)
           / (sqrt(v0 + v1 + v2 + v4 + 2 * sqrt(v0 * v3 + pow(v0 + v1 + v2 + v4, 2) / 4.))
              * sqrt(pow(a, 4) + 2 * v1 * v2 - 2 * (v1 + v2) * v3 + 2 * v0 * (v1 + v2 + v3)
                     + pow(X, 4) + pow(Y, 4) + pow(Z, 4))));
}

inline auto dH_KS_dZ(double M, double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Z, 5)};
  const double v5{pow(Z, 3)};
  const double v6{pow(Y, 4)};
  const double v7{pow(X, 4)};
  const double v8{pow(a, 4)};
  const double v9{pow(Z, 4)};
  const double v10{sqrt(2)};
  const double v11{v0 * v3};
  const double v12{v0 + v1 + v2 - v3};
  const double v13{sqrt(v11 + pow(v12, 2) / 4.)};
  const double v14{sqrt(v12 / 2. + v13)};

  return (-2 * M * pow(v10, 2) * v14
          * (3 * v2 * v4 + v5 * (3 * (v6 + v7) + v8) + 2 * pow(a, 6) * Z
             + (3 * v2 * v6 - 3 * v2 * v8 + pow(X, 6) + pow(Y, 6)) * Z + pow(Z, 7)
             + v1
                   * (3 * v4 + 4 * v13 * v5 + 6 * v2 * v5
                      + (4 * v13 * v2 + 2 * v13 * v3 + 3 * v7 - 3 * v8) * Z)
             + 2 * v13 * (v4 + 2 * v2 * v5 - v3 * v5 + (v2 * v3 + v6 + v7 - 2 * v8) * Z)))
         / (sqrt(2 * v11 + 2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v6 + v7 + v8
                 + v9)
            * pow(2 * v11
                      + 2
                            * (v1 * v13 + v1 * v2 + v13 * v2 + v0 * (v1 + v13 + v2)
                               - (v1 + v13 + v2) * v3)
                      + v6 + v7 + v8 + v9,
                  2));
}

inline auto dl1_KS_dZ(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};
  const double v5{sqrt(2)};
  const double v6{v0 * v3};
  const double v7{v0 + v1 + v2 + v4};
  const double v8{sqrt(4 * v6 + pow(v7, 2))};
  const double v9{sqrt(pow(v0 + v1 + v2 + v4, 2) + 4 * v6) + v7};

  return -(((v5 * X * (v0 + v1 - 3 * v3 + v8 + pow(X, 2)) + 4 * a * sqrt(v9) * Y) * Z)
           / ((v0 + v1 + v2 + v3 + v8) * sqrt(v9)
              * sqrt(pow(a, 4) + 2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + 2 * v6
                     + pow(X, 4) + pow(Y, 4) + pow(Z, 4))));
}

inline auto dl2_KS_dZ(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{-v3};
  const double v5{sqrt(2)};
  const double v6{v0 * v3};
  const double v7{v0 + v1 + v2 + v4};
  const double v8{sqrt(4 * v6 + pow(v7, 2))};
  const double v9{sqrt(pow(v0 + v1 + v2 + v4, 2) + 4 * v6) + v7};

  return ((4 * a * sqrt(v9) * X - v5 * Y * (v0 + v2 - 3 * v3 + v8 + pow(Y, 2))) * Z)
         / ((v0 + v1 + v2 + v3 + v8) * sqrt(v9)
            * sqrt(pow(a, 4) + 2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + 2 * v6
                   + pow(X, 4) + pow(Y, 4) + pow(Z, 4)));
}

inline auto dl3_KS_dZ(double a, double X, double Y, double Z) {
  const double v0{pow(Z, 2)};
  const double v1{pow(Y, 2)};
  const double v2{pow(X, 2)};
  const double v3{pow(a, 2)};
  const double v4{pow(Y, 4)};
  const double v5{pow(X, 4)};
  const double v6{pow(a, 4)};
  const double v7{-v3};
  const double v8{v0 * v3};
  const double v9{v0 + v1 + v2 + v7};
  const double v10{sqrt(v8 + pow(v9, 2) / 4.)};

  return (sqrt(2)
          * (2 * (v1 * v10 + v1 * v2 + v10 * v2 - (v1 + v10 + v2) * v3) + v4 + v5 + v6
             + v0 * (v1 + v2 - v7)))
         / (pow(sqrt(pow(v0 + v1 + v2 + v7, 2) + 4 * v8) + v9, 1.5)
            * sqrt(2 * v1 * v2 + 2 * v0 * (v1 + v2) - 2 * (v1 + v2) * v3 + v4 + v5 + v6 + 2 * v8
                   + pow(Z, 4)));
}

} // namespace sks_aux

#endif // GRLENSING_SKS_AUX_FUNCTIONS