#ifndef GRLENSING_SKS_AUX_FUNCTIONS
#define GRLENSING_SKS_AUX_FUNCTIONS

#include "../../../../include/tensor_types.hpp"

#include <cmath>
#include <cstddef>

namespace sks_aux {

inline auto r_KS(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;

  return sqrt(v0 + v1 + v2 + v4 + sqrt(4 * v0 * v3 + pow(v0 + v1 + v2 + v4, 2))) / sqrt(2);
}

inline auto dr_KS_dX(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * X) / (sqrt(2) * sqrt(v5 + pow(v6, 2)));
}

inline auto dr_KS_dY(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * Y) / (sqrt(2) * sqrt(v5 + pow(v6, 2)));
}

inline auto dr_KS_dZ(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;
  const double v7 = v5 + pow(v6, 2);

  return ((v0 + v1 + v2 + v3 + sqrt(v7)) * Z)
         / (sqrt(2) * sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * sqrt(v7));
}

inline auto H_KS(double M, double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return (sqrt(2) * M * sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6)) / sqrt(v5 + pow(v6, 2));
}

inline auto dH_KS_dX(double M, double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;
  const double v7 = v5 + pow(v6, 2);
  const double v8 = sqrt(v7);

  return -((sqrt(2) * M * X
            * (pow(a, 4) + 2 * v0 * v1 + 2 * v0 * v2 + 2 * v1 * v2 - 6 * v0 * v3 - 2 * v1 * v3
               - 2 * v2 * v3 + v4 * sqrt(pow(v0 + v1 + v2 - v3, 2) + v5) + v0 * v8 + v1 * v8
               + v2 * v8 + pow(X, 4) + pow(Y, 4) + pow(Z, 4)))
           / (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * pow(v7, 1.5)));
}

inline auto dH_KS_dY(double M, double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;
  const double v7 = v5 + pow(v6, 2);
  const double v8 = sqrt(v7);

  return -((sqrt(2) * M * Y
            * (pow(a, 4) + 2 * v0 * v1 + 2 * v0 * v2 + 2 * v1 * v2 - 6 * v0 * v3 - 2 * v1 * v3
               - 2 * v2 * v3 + v4 * sqrt(pow(v0 + v1 + v2 - v3, 2) + v5) + v0 * v8 + v1 * v8
               + v2 * v8 + pow(X, 4) + pow(Y, 4) + pow(Z, 4)))
           / (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * pow(v7, 1.5)));
}

inline auto dH_KS_dZ(double M, double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;
  const double v7 = v5 + pow(v6, 2);
  const double v8 = sqrt(v7);

  return -((sqrt(2) * M * Z
            * (-3 * pow(a, 4) + 2 * v0 * v1 + 2 * v0 * v2 + 2 * v1 * v2 - 2 * v0 * v3 + 2 * v1 * v3
               + 2 * v2 * v3 + v0 * v8 + v1 * v8 + v2 * v8 + v3 * v8 + pow(X, 4) + pow(Y, 4)
               + pow(Z, 4)))
           / (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * pow(v7, 1.5)));
}

inline auto l1_KS(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return (sqrt(2) * sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * X + 2 * a * Y)
         / (v0 + v1 + v2 + v3 + sqrt(v5 + pow(v6, 2)));
}

inline auto dl1_KS_dX(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = sqrt(2);
  const double v6 = 4 * v0 * v3;
  const double v7 = v0 + v1 + v2 + v4;
  const double v8 = v6 + pow(v7, 2);
  const double v9 = sqrt(v8);
  const double v10 = sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v6) + v7);

  return (pow(a, 4) * v10 * v5 + 2 * v0 * v1 * v10 * v5 + v0 * v10 * v2 * v5 + v1 * v10 * v2 * v5
          + 2 * v0 * v10 * v3 * v5 - 2 * v1 * v10 * v3 * v5 + v10 * v2 * v3 * v5
          + v0 * v10 * v5 * v9 + v1 * v10 * v5 * v9 + v10 * v3 * v5 * v9 + 4 * pow(a, 3) * X * Y
          - 4 * a * v0 * X * Y - 4 * a * v9 * X * Y - 4 * a * pow(X, 3) * Y - 4 * a * X * pow(Y, 3)
          + v10 * v5 * pow(Y, 4) + v10 * v5 * pow(Z, 4))
         / (sqrt(v8) * pow(v0 + v1 + v2 + v3 + v9, 2));
}

inline auto dl1_KS_dY(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = pow(a, 3);
  const double v5 = -v3;
  const double v6 = sqrt(2);
  const double v7 = 4 * v0 * v3;
  const double v8 = v0 + v1 + v2 + v5;
  const double v9 = v7 + pow(v8, 2);
  const double v10 = sqrt(v9);
  const double v11 = sqrt(sqrt(pow(v0 + v1 + v2 + v5, 2) + v7) + v8);

  return (2 * pow(a, 5) + 2 * a * v0 * v10 - 2 * a * v1 * v10 + 4 * a * v0 * v2 + 2 * a * v10 * v2
          + 4 * v0 * v4 + 2 * v10 * v4 - 4 * v2 * v4 + 2 * a * pow(X, 4) - v0 * v11 * v6 * X * Y
          - v10 * v11 * v6 * X * Y + 3 * v11 * v3 * v6 * X * Y - v11 * v6 * pow(X, 3) * Y
          - v11 * v6 * X * pow(Y, 3) - 2 * a * pow(Y, 4) + 2 * a * pow(Z, 4))
         / (pow(v0 + v1 + v10 + v2 + v3, 2) * sqrt(v9));
}

inline auto dl1_KS_dZ(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = sqrt(2);
  const double v6 = 4 * v0 * v3;
  const double v7 = v0 + v1 + v2 + v4;
  const double v8 = v6 + pow(v7, 2);
  const double v9 = sqrt(v8);
  const double v10 = sqrt(pow(v0 + v1 + v2 + v4, 2) + v6) + v7;

  return -(((v0 * v5 * X + v1 * v5 * X - 3 * v3 * v5 * X + v5 * v9 * X + v5 * pow(X, 3)
             + 4 * a * sqrt(v10) * Y)
            * Z)
           / (sqrt(v10) * sqrt(v8) * (v0 + v1 + v2 + v3 + v9)));
}

inline auto l2_KS(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return (-2 * a * X + sqrt(2) * sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * Y)
         / (v0 + v1 + v2 + v3 + sqrt(v5 + pow(v6, 2)));
}

inline auto dl2_KS_dX(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = pow(a, 3);
  const double v5 = -v3;
  const double v6 = sqrt(2);
  const double v7 = 4 * v0 * v3;
  const double v8 = v0 + v1 + v2 + v5;
  const double v9 = v7 + pow(v8, 2);
  const double v10 = sqrt(v9);
  const double v11 = sqrt(sqrt(pow(v0 + v1 + v2 + v5, 2) + v7) + v8);

  return (-2 * pow(a, 5) - 4 * a * v0 * v1 - 2 * a * v0 * v10 - 2 * a * v1 * v10 + 2 * a * v10 * v2
          - 4 * v0 * v4 + 4 * v1 * v4 - 2 * v10 * v4 + 2 * a * pow(X, 4) - v0 * v11 * v6 * X * Y
          - v10 * v11 * v6 * X * Y + 3 * v11 * v3 * v6 * X * Y - v11 * v6 * pow(X, 3) * Y
          - v11 * v6 * X * pow(Y, 3) - 2 * a * pow(Y, 4) - 2 * a * pow(Z, 4))
         / (pow(v0 + v1 + v10 + v2 + v3, 2) * sqrt(v9));
}

inline auto dl2_KS_dY(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = sqrt(2);
  const double v6 = 4 * v0 * v3;
  const double v7 = v0 + v1 + v2 + v4;
  const double v8 = v6 + pow(v7, 2);
  const double v9 = sqrt(v8);
  const double v10 = sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v6) + v7);

  return (pow(a, 4) * v10 * v5 + v0 * v1 * v10 * v5 + 2 * v0 * v10 * v2 * v5 + v1 * v10 * v2 * v5
          + 2 * v0 * v10 * v3 * v5 + v1 * v10 * v3 * v5 - 2 * v10 * v2 * v3 * v5
          + v0 * v10 * v5 * v9 + v10 * v2 * v5 * v9 + v10 * v3 * v5 * v9 + v10 * v5 * pow(X, 4)
          - 4 * pow(a, 3) * X * Y + 4 * a * v0 * X * Y + 4 * a * v9 * X * Y + 4 * a * pow(X, 3) * Y
          + 4 * a * X * pow(Y, 3) + v10 * v5 * pow(Z, 4))
         / (sqrt(v8) * pow(v0 + v1 + v2 + v3 + v9, 2));
}

inline auto dl2_KS_dZ(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = sqrt(2);
  const double v6 = 4 * v0 * v3;
  const double v7 = v0 + v1 + v2 + v4;
  const double v8 = v6 + pow(v7, 2);
  const double v9 = sqrt(v8);
  const double v10 = sqrt(pow(v0 + v1 + v2 + v4, 2) + v6) + v7;

  return -(((-4 * a * sqrt(v10) * X + v0 * v5 * Y + v2 * v5 * Y - 3 * v3 * v5 * Y + v5 * v9 * Y
             + v5 * pow(Y, 3))
            * Z)
           / (sqrt(v10) * sqrt(v8) * (v0 + v1 + v2 + v3 + v9)));
}

inline auto l3_KS(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;

  return (sqrt(2) * Z) / sqrt(v0 + v1 + v2 + v4 + sqrt(4 * v0 * v3 + pow(v0 + v1 + v2 + v4, 2)));
}

inline auto dl3_KS_dX(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return -((sqrt(2) * X * Z)
           / (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * sqrt(v5 + pow(v6, 2))));
}

inline auto dl3_KS_dY(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;

  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;

  return -((sqrt(2) * Y * Z)
           / (sqrt(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6) * sqrt(v5 + pow(v6, 2))));
}

inline auto dl3_KS_dZ(double a, double X, double Y, double Z) -> double {
  using std::pow;
  using std::sqrt;
  const double v0 = pow(Z, 2);
  const double v1 = pow(Y, 2);
  const double v2 = pow(X, 2);
  const double v3 = pow(a, 2);
  const double v4 = -v3;
  const double v5 = 4 * v0 * v3;
  const double v6 = v0 + v1 + v2 + v4;
  const double v7 = v5 + pow(v6, 2);
  const double v8 = sqrt(v7);

  return -((sqrt(2)
            * (-pow(a, 4) - v0 * v1 - v0 * v2 - 2 * v1 * v2 + 2 * v1 * v3 + 2 * v2 * v3 + v0 * v4
               - v1 * v8 - v2 * v8 + v3 * v8 - pow(X, 4) - pow(Y, 4)))
           / (pow(sqrt(pow(v0 + v1 + v2 + v4, 2) + v5) + v6, 1.5) * sqrt(v7)));
}

} // namespace sks_aux

#endif // GRLENSING_SKS_AUX_FUNCTIONS