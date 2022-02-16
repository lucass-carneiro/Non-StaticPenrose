
#include "KerrSchild_Kerr.hpp"

using grlensing::kernel;
using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;

GRLENSING_KERRSCHILD_KERR_METRIC_API void
KerrSchild_Kerr::load_parameters(const YAML::Node &) {}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::lapse(double, double,
                                                                 double, double)
    -> double {
  return 1.0;
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::l_shift(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{0.0, 0.0, 0.0};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::u_shift(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{0.0, 0.0, 0.0};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::ll_smetric(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{
      {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::uu_smetric(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{
      {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::ll_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::ul_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::spatial_christoffel(double, double, double, double)
    -> metric_server::chirstofell_t {
  metric_server::chirstofell_t Gamma;

  Gamma[0][0][0] = 0.0;
  Gamma[0][0][1] = 0.0;
  Gamma[0][0][2] = 0.0;
  Gamma[0][1][1] = 0.0;
  Gamma[0][1][2] = 0.0;
  Gamma[0][2][2] = 0.0;

  Gamma[1][0][0] = 0.0;
  Gamma[1][0][1] = 0.0;
  Gamma[1][0][2] = 0.0;
  Gamma[1][1][1] = 0.0;
  Gamma[1][1][2] = 0.0;
  Gamma[1][2][2] = 0.0;

  Gamma[2][0][0] = 0.0;
  Gamma[2][0][1] = 0.0;
  Gamma[2][0][2] = 0.0;
  Gamma[2][1][1] = 0.0;
  Gamma[2][1][2] = 0.0;
  Gamma[2][2][2] = 0.0;

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

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::grad_lapse(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{0.0, 0.0, 0.0};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto
KerrSchild_Kerr::grad_ushift(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{
      {{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::name()
    -> std::string_view {
  return "Kerr-Schild Kerr";
}

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API auto get_engine_version()
    -> unsigned {
  return unsigned(1);
}

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API void
register_plugin(kernel &kernel) {
  kernel.get_metric_server().add_metric(
      std::unique_ptr<metric_server::adm_metric>(new KerrSchild_Kerr()));
}
