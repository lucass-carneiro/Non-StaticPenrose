
#include "SKS.hpp"

using grlensing::kernel;
using grlensing::metric_server;
using grlensing::SKS;

GRLENSING_SKS_METRIC_API void SKS::load_parameters(const YAML::Node &node) {

  M1 = node["SKS_Settings"]["M1"].as<double>();
  a1 = node["SKS_Settings"]["a1"].as<double>();

  M2 = node["SKS_Settings"]["M2"].as<double>();
  a2 = node["SKS_Settings"]["a2"].as<double>();
}

GRLENSING_SKS_METRIC_API auto SKS::lapse(double, double, double, double) -> double { return 1.0; }

GRLENSING_SKS_METRIC_API auto SKS::l_shift(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{};
}

GRLENSING_SKS_METRIC_API auto SKS::u_shift(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{};
}

GRLENSING_SKS_METRIC_API auto SKS::ll_smetric(double, double, double, double)
    -> metric_server::spatial_matrix {
  metric_server::spatial_matrix metric{};

  metric[0][0] = 1.0;
  metric[1][1] = metric[0][0];
  metric[2][2] = metric[0][0];

  return metric;
}

GRLENSING_SKS_METRIC_API auto SKS::uu_smetric(double, double, double, double)
    -> metric_server::spatial_matrix {
  metric_server::spatial_matrix metric{};

  metric[0][0] = 1.0;
  metric[1][1] = metric[0][0];
  metric[2][2] = metric[0][0];

  return metric;
}

GRLENSING_SKS_METRIC_API auto SKS::ll_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{};
}

GRLENSING_SKS_METRIC_API auto SKS::ul_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{};
}

GRLENSING_SKS_METRIC_API auto SKS::spatial_christoffel(double, double, double, double)
    -> metric_server::chirstofell_t {
  return metric_server::chirstofell_t{};
}

GRLENSING_SKS_METRIC_API auto SKS::grad_lapse(double, double, double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{};
}

GRLENSING_SKS_METRIC_API auto SKS::grad_ushift(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{};
}

GRLENSING_SKS_METRIC_API auto SKS::name() -> std::string_view {
  return "Superimposed Kerr binary in Kerr-Schild coordinates";
}

extern "C" GRLENSING_SKS_METRIC_API auto get_engine_version() -> unsigned { return unsigned(1); }

extern "C" GRLENSING_SKS_METRIC_API void register_plugin(kernel &kernel) {
  kernel.get_metric_server().add_metric(std::unique_ptr<metric_server::adm_metric>(new SKS()));
}
