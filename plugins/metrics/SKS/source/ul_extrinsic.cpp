#include "SKS.hpp"
#include "aux_functions.hpp"

#include <cstddef>

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::ul_extrinsic(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {
  using namespace sks_aux;

  const auto llk = ll_extrinsic(t, x, y, z);
  const auto uugamma = uu_smetric(t, x, y, z);

  metric_server::spatial_matrix ulK{};

  for (size_t i = 0; i <= 2; i++) {
    for (size_t j = 0; j <= 2; j++) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
      ulK[i][j] = uugamma(i, 0) * llk(j, 0) + uugamma(i, 1) * llk(j, 1) + uugamma(i, 2) * llk(j, 2);
    }
  }

  return ulK;
}