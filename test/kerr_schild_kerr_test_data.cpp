#include "kerr_schild_kerr_test_data.hpp"

#include "test_constants.hpp"

auto grlensing_tests::r(double a, double x, double y, double z) -> double {
  double part1 = Power(x, 2) + Power(y, 2) + Power(z, 2) - Power(a, 2);
  return Sqrt(part1 + Sqrt(4 * Power(a, 2) * Power(z, 2) + Power(part1, 2))) / Sqrt(2);
}

auto grlensing_tests::H(double M, double a, double x, double y, double z) -> double {
  double rVal = r(a, x, y, z);
  return (M * rVal) / (Power(rVal, 2) + Power(a, 2) * Power(z / rVal, 2));
}

auto grlensing_tests::l(double a, double x, double y, double z)
    -> grlensing::callable_array<4, double> {

  grlensing::callable_array<4, double> kvec{};

  double rVal = r(a, x, y, z);
  double den = Power(rVal, 2) + Power(a, 2);

  kvec[0] = 1.0;
  kvec[1] = (rVal * x + a * y) / den;
  kvec[2] = (rVal * y - a * x) / den;
  kvec[3] = z / rVal;

  return kvec;
}

auto grlensing_tests::ll_g(double M, double a, double x, double y, double z)
    -> grlensing::callable_matrix<4, 4, double> {

  auto Hval = H(M, a, x, y, z);
  auto lVal = l(a, x, y, z);

  grlensing::callable_matrix<4, 4, double> llg{};

  llg[0][0] = 2 * Hval * lVal[0] * lVal[0] - 1.0;
  llg[0][1] = 2 * Hval * lVal[0] * lVal[1];
  llg[0][2] = 2 * Hval * lVal[0] * lVal[2];
  llg[0][3] = 2 * Hval * lVal[0] * lVal[3];

  llg[1][1] = 2 * Hval * lVal[1] * lVal[1] + 1.0;
  llg[1][2] = 2 * Hval * lVal[1] * lVal[2];
  llg[1][3] = 2 * Hval * lVal[1] * lVal[3];

  llg[2][2] = 2 * Hval * lVal[2] * lVal[2] + 1.0;
  llg[2][3] = 2 * Hval * lVal[2] * lVal[3];

  llg[3][3] = 2 * Hval * lVal[3] * lVal[3] + 1.0;

  llg[1][0] = llg[0][1];
  llg[2][0] = llg[0][2];
  llg[2][1] = llg[1][2];
  llg[3][0] = llg[0][3];
  llg[3][1] = llg[1][3];
  llg[3][2] = llg[2][3];

  return llg;
}

auto grlensing_tests::uu_g(double M, double a, double x, double y, double z)
    -> grlensing::callable_matrix<4, 4, double> {

  const auto llg = ll_g(M, a, x, y, z);
  grlensing::callable_matrix<4, 4, double> uug{};

  uug[0][0]
      = (Power(llg(1, 3), 2) * llg(2, 2) - 2 * llg(1, 2) * llg(1, 3) * llg(2, 3)
         + Power(llg(1, 2), 2) * llg(3, 3)
         + llg(1, 1) * (Power(llg(2, 3), 2) - llg(2, 2) * llg(3, 3)))
        / (llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2)
           + Power(llg(0, 3), 2) * (-Power(llg(1, 2), 2) + llg(1, 1) * llg(2, 2))
           - 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           - Power(llg(0, 1), 2) * Power(llg(2, 3), 2) + llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (llg(1, 2) * llg(1, 3) - llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
           + (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (-Power(llg(1, 3), 2) + llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)));
  uug[0][1]
      = (llg(0, 3) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3))
         + llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3))
         + llg(0, 1) * (-Power(llg(2, 3), 2) + llg(2, 2) * llg(3, 3)))
        / (llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2)
           + Power(llg(0, 3), 2) * (-Power(llg(1, 2), 2) + llg(1, 1) * llg(2, 2))
           - 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           - Power(llg(0, 1), 2) * Power(llg(2, 3), 2) + llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (llg(1, 2) * llg(1, 3) - llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
           + (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (-Power(llg(1, 3), 2) + llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)));
  uug[0][2]
      = (llg(0, 3) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
         + llg(0, 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
         + llg(0, 1) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)))
        / (-(llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2))
           + Power(llg(0, 3), 2) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
           + 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           + Power(llg(0, 1), 2) * Power(llg(2, 3), 2) - llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
           - (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)));
  uug[0][3]
      = (llg(0, 3) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
         + llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
         + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
        / (-(llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2))
           + Power(llg(0, 3), 2) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
           + 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           + Power(llg(0, 1), 2) * Power(llg(2, 3), 2) - llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
           - (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)));

  uug[1][1]
      = (Power(llg(0, 3), 2) * llg(2, 2) - 2 * llg(0, 2) * llg(0, 3) * llg(2, 3)
         + Power(llg(0, 2), 2) * llg(3, 3)
         + llg(0, 0) * (Power(llg(2, 3), 2) - llg(2, 2) * llg(3, 3)))
        / (llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2)
           + Power(llg(0, 3), 2) * (-Power(llg(1, 2), 2) + llg(1, 1) * llg(2, 2))
           - 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           - Power(llg(0, 1), 2) * Power(llg(2, 3), 2) + llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (llg(1, 2) * llg(1, 3) - llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
           + (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (-Power(llg(1, 3), 2) + llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)));
  uug[1][2]
      = (Power(llg(0, 3), 2) * llg(1, 2)
         - llg(0, 3) * (llg(0, 2) * llg(1, 3) + llg(0, 1) * llg(2, 3))
         + llg(0, 1) * llg(0, 2) * llg(3, 3)
         + llg(0, 0) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)))
        / (-(llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2))
           + Power(llg(0, 3), 2) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
           + 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           + Power(llg(0, 1), 2) * Power(llg(2, 3), 2) - llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
           - (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)));
  uug[1][3]
      = (Power(llg(0, 2), 2) * llg(1, 3) + llg(0, 1) * llg(0, 3) * llg(2, 2)
         - llg(0, 2) * (llg(0, 3) * llg(1, 2) + llg(0, 1) * llg(2, 3))
         + llg(0, 0) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
        / (-(llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2))
           + Power(llg(0, 3), 2) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
           + 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           + Power(llg(0, 1), 2) * Power(llg(2, 3), 2) - llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
           - (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)));

  uug[2][2]
      = (Power(llg(0, 3), 2) * llg(1, 1) - 2 * llg(0, 1) * llg(0, 3) * llg(1, 3)
         + Power(llg(0, 1), 2) * llg(3, 3)
         + llg(0, 0) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3)))
        / (llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2)
           + Power(llg(0, 3), 2) * (-Power(llg(1, 2), 2) + llg(1, 1) * llg(2, 2))
           - 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           - Power(llg(0, 1), 2) * Power(llg(2, 3), 2) + llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (llg(1, 2) * llg(1, 3) - llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
           + (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (-Power(llg(1, 3), 2) + llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)));
  uug[2][3]
      = (llg(1, 2) * (-(llg(0, 1) * llg(0, 3)) + llg(0, 0) * llg(1, 3))
         + llg(0, 2) * (llg(0, 3) * llg(1, 1) - llg(0, 1) * llg(1, 3))
         + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 3))
        / (-(llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2))
           + Power(llg(0, 3), 2) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2))
           + 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           + Power(llg(0, 1), 2) * Power(llg(2, 3), 2) - llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (-(llg(1, 2) * llg(1, 3)) + llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (llg(1, 3) * llg(2, 2) - llg(1, 2) * llg(2, 3)))
           - (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (Power(llg(1, 3), 2) - llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (-(llg(1, 3) * llg(2, 3)) + llg(1, 2) * llg(3, 3)));

  uug[3][3]
      = (Power(llg(0, 2), 2) * llg(1, 1) - 2 * llg(0, 1) * llg(0, 2) * llg(1, 2)
         + Power(llg(0, 1), 2) * llg(2, 2)
         + llg(0, 0) * (Power(llg(1, 2), 2) - llg(1, 1) * llg(2, 2)))
        / (llg(0, 0) * Power(llg(1, 3), 2) * llg(2, 2)
           + Power(llg(0, 3), 2) * (-Power(llg(1, 2), 2) + llg(1, 1) * llg(2, 2))
           - 2 * llg(0, 0) * llg(1, 2) * llg(1, 3) * llg(2, 3)
           - Power(llg(0, 1), 2) * Power(llg(2, 3), 2) + llg(0, 0) * llg(1, 1) * Power(llg(2, 3), 2)
           + 2 * llg(0, 3)
                 * (llg(0, 2) * (llg(1, 2) * llg(1, 3) - llg(1, 1) * llg(2, 3))
                    + llg(0, 1) * (-(llg(1, 3) * llg(2, 2)) + llg(1, 2) * llg(2, 3)))
           + (llg(0, 0) * Power(llg(1, 2), 2)
              + (Power(llg(0, 1), 2) - llg(0, 0) * llg(1, 1)) * llg(2, 2))
                 * llg(3, 3)
           + Power(llg(0, 2), 2) * (-Power(llg(1, 3), 2) + llg(1, 1) * llg(3, 3))
           + 2 * llg(0, 1) * llg(0, 2) * (llg(1, 3) * llg(2, 3) - llg(1, 2) * llg(3, 3)));

  uug[1][0] = uug[0][1];
  uug[2][0] = uug[0][2];
  uug[2][1] = uug[1][2];
  uug[3][0] = uug[0][3];
  uug[3][1] = uug[1][3];
  uug[3][2] = uug[2][3];

  return uug;
}