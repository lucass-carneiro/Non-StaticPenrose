
#include "SKS.hpp"

using grlensing::kernel;
using grlensing::metric_server;
using grlensing::SKS;

GRLENSING_SKS_METRIC_API void SKS::load_parameters(const YAML::Node &node) {

  M1 = node["SKS_Settings"]["M1"].as<double>();
  a1 = node["SKS_Settings"]["a1"].as<double>();

  M2 = node["SKS_Settings"]["M2"].as<double>();
  a2 = node["SKS_Settings"]["a2"].as<double>();

  b = node["SKS_Settings"]["b"].as<double>();
}

GRLENSING_SKS_METRIC_API auto SKS::name() -> std::string_view {
  return "Superimposed Kerr binary in Kerr-Schild coordinates";
}

extern "C" GRLENSING_SKS_METRIC_API auto get_engine_version() -> unsigned { return unsigned(1); }

extern "C" GRLENSING_SKS_METRIC_API void register_plugin(kernel &kernel) {
  kernel.get_metric_server().add_metric(std::unique_ptr<metric_server::adm_metric>(new SKS()));
}
