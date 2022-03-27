#include "KerrSchild_Kerr.hpp"

#include "aux_functions.hpp"

using grlensing::kernel;
using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;

GRLENSING_KERRSCHILD_KERR_METRIC_API void KerrSchild_Kerr::load_parameters(const YAML::Node &node) {

  M = node["KerrSchild_Kerr_Settings"]["M"].as<double>();
  a = node["KerrSchild_Kerr_Settings"]["a"].as<double>();

  double delta = M * M - a * a;

  if (delta < 0.0) {
    log<LogEvent::warning>("The spin parameter {} together with the mass {} produces a naked "
                           "singularity. Some quantities might be ill defined.",
                           a, M);
  }
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::name() -> std::string_view {
  return "Kerr-Schild Kerr";
}

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API auto get_engine_version() -> unsigned {
  return unsigned(1);
}

extern "C" GRLENSING_KERRSCHILD_KERR_METRIC_API void register_plugin(kernel &kernel) {
  kernel.get_metric_server().add_metric(
      std::unique_ptr<metric_server::adm_metric>(new KerrSchild_Kerr()));
}