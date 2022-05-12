#include "SKS.hpp"
#include "aux_functions.hpp"

#include <cstddef>

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::ll_extrinsic(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {
  using namespace sks_aux;

  const auto alpha = lapse(t, x, y, z);
  const auto Gamma = spatial_christoffel(t, x, y, z);
  const auto lbeta = l_shift(t, x, y, z);
  const auto dllg = dllg_SKS(M1, M2, a1, a2, b, t, x, y, z);

  metric_server::spatial_matrix llK{};

  for (size_t i = 0; i <= 2; i++) {
    for (size_t j = 0; j <= 2; j++) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
      llK[i][j] = (-dllg(0, i + 1, j + 1)
                   + (dllg(i + 1, 0, j + 1) - Gamma(0, i, j) * lbeta(0) - Gamma(1, i, j) * lbeta(1)
                      - Gamma(2, i, j) * lbeta(2))
                   + (dllg(j + 1, 0, i + 1) - Gamma(0, j, i) * lbeta(0) - Gamma(1, j, i) * lbeta(1)
                      - Gamma(2, j, i) * lbeta(2)))
                  / (2 * alpha);
    }
  }

  return llK;
}