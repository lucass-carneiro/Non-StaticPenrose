
#ifndef GRLENSING_KERRSCHILD_KERR_METRIC_HPP
#define GRLENSING_KERRSCHILD_KERR_METRIC_HPP

#include "api_macros.hpp"

namespace grlensing {

class KerrSchild_Kerr final : public metric_server::adm_metric {
public:
  GRLENSING_KERRSCHILD_KERR_METRIC_API void load_parameters(const YAML::Node &) final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto lapse(double, double, double, double) -> double final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto l_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto u_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto ll_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto uu_smetric(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto ll_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto ul_extrinsic(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto spatial_christoffel(double, double, double, double)
      -> metric_server::chirstofell_t final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto grad_lapse(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto grad_ushift(double, double, double, double)
      -> metric_server::spatial_matrix final;

  GRLENSING_KERRSCHILD_KERR_METRIC_API auto name() -> std::string_view final;

private:
  double M = 1.0;
  double a = M / 2;
};

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API auto get_engine_version() -> unsigned;

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API void register_plugin(kernel &kernel);

} // namespace grlensing

#endif // GRLENSING_KERRSCHILD_KERR_METRIC_HPP ";
