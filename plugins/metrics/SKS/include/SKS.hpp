#ifndef GRLENSING_SKS_METRIC_HPP
#define GRLENSING_SKS_METRIC_HPP

#include "api_macros.hpp"

namespace grlensing {

class SKS final : public metric_server::adm_metric {
public:
  GRLENSING_SKS_METRIC_API void load_parameters(const YAML::Node &) final;

  GRLENSING_SKS_METRIC_API auto lapse(double, double, double, double) -> double final;

  GRLENSING_SKS_METRIC_API auto l_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto u_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto ll_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto uu_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto ll_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto ul_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto spatial_christoffel(double, double, double, double)
      -> metric_server::chirstofell_t final;

  GRLENSING_SKS_METRIC_API auto grad_lapse(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_SKS_METRIC_API auto grad_ushift(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_SKS_METRIC_API auto name() -> std::string_view final;

private:
  double M1 = 0.5;
  double M2 = 0.5;
  double a1 = 0.9 * M1;
  double a2 = 0.9 * M2;
  double b = 20.0;
};

extern "C" GRLENSING_SKS_METRIC_API auto get_engine_version() -> unsigned;

extern "C" GRLENSING_SKS_METRIC_API void register_plugin(kernel &kernel);

} // namespace grlensing

#endif // GRLENSING_SKS_METRIC_HPP ";
