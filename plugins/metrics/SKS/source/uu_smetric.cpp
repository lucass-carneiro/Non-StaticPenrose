#include "SKS.hpp"
#include "aux_functions.hpp"

using namespace grlensing;

GRLENSING_SKS_METRIC_API auto SKS::uu_smetric(double t, double x, double y, double z)
    -> metric_server::spatial_matrix {
  using namespace sks_aux;

  const auto llgamma = ll_smetric(t, x, y, z);
  const auto uugamma = inverse_symmetric_3x3(llgamma);

  return uugamma;
}