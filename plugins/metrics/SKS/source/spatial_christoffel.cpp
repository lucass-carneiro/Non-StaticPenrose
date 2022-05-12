#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::spatial_christoffel(double t, double x, double y, double z)
    -> metric_server::chirstofell_t {

  using namespace sks_aux;

  const auto uugamma = uu_smetric(t, x, y, z);
  const auto dllg = dllg_SKS(M1, M2, a1, a2, b, t, x, y, z);

  metric_server::chirstofell_t Gamma{};

  Gamma[0][0][0]
      = (dllg(1, 1, 1) * uugamma(0, 0) + (2 * dllg(1, 2, 1) - dllg(2, 1, 1)) * uugamma(0, 1)
         + (2 * dllg(1, 3, 1) - dllg(3, 1, 1)) * uugamma(0, 2))
        / 2.;
  Gamma[0][0][1] = (dllg(2, 1, 1) * uugamma(0, 0)
                    + (dllg(1, 2, 2) - dllg(2, 1, 2) + dllg(2, 2, 1)) * uugamma(0, 1)
                    + (dllg(1, 3, 2) + dllg(2, 3, 1) - dllg(3, 1, 2)) * uugamma(0, 2))
                   / 2.;
  Gamma[0][0][2] = (dllg(3, 1, 1) * uugamma(0, 0)
                    + (dllg(1, 2, 3) - dllg(2, 1, 3) + dllg(3, 2, 1)) * uugamma(0, 1)
                    + (dllg(1, 3, 3) - dllg(3, 1, 3) + dllg(3, 3, 1)) * uugamma(0, 2))
                   / 2.;
  Gamma[0][1][1]
      = (-((dllg(1, 2, 2) - 2 * dllg(2, 1, 2)) * uugamma(0, 0)) + dllg(2, 2, 2) * uugamma(0, 1)
         + (2 * dllg(2, 3, 2) - dllg(3, 2, 2)) * uugamma(0, 2))
        / 2.;
  Gamma[0][1][2] = ((-dllg(1, 2, 3) + dllg(2, 1, 3) + dllg(3, 1, 2)) * uugamma(0, 0)
                    + dllg(3, 2, 2) * uugamma(0, 1)
                    + (dllg(2, 3, 3) - dllg(3, 2, 3) + dllg(3, 3, 2)) * uugamma(0, 2))
                   / 2.;
  Gamma[0][2][2]
      = (-((dllg(1, 3, 3) - 2 * dllg(3, 1, 3)) * uugamma(0, 0))
         - (dllg(2, 3, 3) - 2 * dllg(3, 2, 3)) * uugamma(0, 1) + dllg(3, 3, 3) * uugamma(0, 2))
        / 2.;

  Gamma[1][0][0]
      = (dllg(1, 1, 1) * uugamma(1, 0) + (2 * dllg(1, 2, 1) - dllg(2, 1, 1)) * uugamma(1, 1)
         + (2 * dllg(1, 3, 1) - dllg(3, 1, 1)) * uugamma(1, 2))
        / 2.;
  Gamma[1][0][1] = (dllg(2, 1, 1) * uugamma(1, 0)
                    + (dllg(1, 2, 2) - dllg(2, 1, 2) + dllg(2, 2, 1)) * uugamma(1, 1)
                    + (dllg(1, 3, 2) + dllg(2, 3, 1) - dllg(3, 1, 2)) * uugamma(1, 2))
                   / 2.;
  Gamma[1][0][2] = (dllg(3, 1, 1) * uugamma(1, 0)
                    + (dllg(1, 2, 3) - dllg(2, 1, 3) + dllg(3, 2, 1)) * uugamma(1, 1)
                    + (dllg(1, 3, 3) - dllg(3, 1, 3) + dllg(3, 3, 1)) * uugamma(1, 2))
                   / 2.;
  Gamma[1][1][1]
      = (-((dllg(1, 2, 2) - 2 * dllg(2, 1, 2)) * uugamma(1, 0)) + dllg(2, 2, 2) * uugamma(1, 1)
         + (2 * dllg(2, 3, 2) - dllg(3, 2, 2)) * uugamma(1, 2))
        / 2.;
  Gamma[1][1][2] = ((-dllg(1, 2, 3) + dllg(2, 1, 3) + dllg(3, 1, 2)) * uugamma(1, 0)
                    + dllg(3, 2, 2) * uugamma(1, 1)
                    + (dllg(2, 3, 3) - dllg(3, 2, 3) + dllg(3, 3, 2)) * uugamma(1, 2))
                   / 2.;
  Gamma[1][2][2]
      = (-((dllg(1, 3, 3) - 2 * dllg(3, 1, 3)) * uugamma(1, 0))
         - (dllg(2, 3, 3) - 2 * dllg(3, 2, 3)) * uugamma(1, 1) + dllg(3, 3, 3) * uugamma(1, 2))
        / 2.;

  Gamma[2][0][0]
      = (dllg(1, 1, 1) * uugamma(2, 0) + (2 * dllg(1, 2, 1) - dllg(2, 1, 1)) * uugamma(2, 1)
         + (2 * dllg(1, 3, 1) - dllg(3, 1, 1)) * uugamma(2, 2))
        / 2.;
  Gamma[2][0][1] = (dllg(2, 1, 1) * uugamma(2, 0)
                    + (dllg(1, 2, 2) - dllg(2, 1, 2) + dllg(2, 2, 1)) * uugamma(2, 1)
                    + (dllg(1, 3, 2) + dllg(2, 3, 1) - dllg(3, 1, 2)) * uugamma(2, 2))
                   / 2.;
  Gamma[2][0][2] = (dllg(3, 1, 1) * uugamma(2, 0)
                    + (dllg(1, 2, 3) - dllg(2, 1, 3) + dllg(3, 2, 1)) * uugamma(2, 1)
                    + (dllg(1, 3, 3) - dllg(3, 1, 3) + dllg(3, 3, 1)) * uugamma(2, 2))
                   / 2.;
  Gamma[2][1][1]
      = (-((dllg(1, 2, 2) - 2 * dllg(2, 1, 2)) * uugamma(2, 0)) + dllg(2, 2, 2) * uugamma(2, 1)
         + (2 * dllg(2, 3, 2) - dllg(3, 2, 2)) * uugamma(2, 2))
        / 2.;
  Gamma[2][1][2] = ((-dllg(1, 2, 3) + dllg(2, 1, 3) + dllg(3, 1, 2)) * uugamma(2, 0)
                    + dllg(3, 2, 2) * uugamma(2, 1)
                    + (dllg(2, 3, 3) - dllg(3, 2, 3) + dllg(3, 3, 2)) * uugamma(2, 2))
                   / 2.;
  Gamma[2][2][2]
      = (-((dllg(1, 3, 3) - 2 * dllg(3, 1, 3)) * uugamma(2, 0))
         - (dllg(2, 3, 3) - 2 * dllg(3, 2, 3)) * uugamma(2, 1) + dllg(3, 3, 3) * uugamma(2, 2))
        / 2.;

  Gamma[0][1][0] = Gamma[0][0][1];
  Gamma[0][2][0] = Gamma[0][0][2];
  Gamma[0][2][1] = Gamma[0][1][2];

  Gamma[1][1][0] = Gamma[1][0][1];
  Gamma[1][2][0] = Gamma[1][0][2];
  Gamma[1][2][1] = Gamma[1][1][2];

  Gamma[2][1][0] = Gamma[2][0][1];
  Gamma[2][2][0] = Gamma[2][0][2];
  Gamma[2][2][1] = Gamma[2][1][2];

  return Gamma;
}