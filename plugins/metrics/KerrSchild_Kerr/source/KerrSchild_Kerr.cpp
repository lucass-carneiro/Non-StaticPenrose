#include "KerrSchild_Kerr.hpp"

#include "aux_functions.hpp"

using grlensing::kernel;
using grlensing::KerrSchild_Kerr;
using grlensing::metric_server;
using ksk_aux::Power;
using ksk_aux::r_KS;
using ksk_aux::Sqrt;

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

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::lapse(double, double x, double y,
                                                                 double z) -> double {

  auto r = r_KS(a, x, y, z);
  return Sqrt(1
              - (2 * M * Power(r, 3) * (Power(a, 2) + Power(r, 2)))
                    / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                       + Power(a, 2) * Power(z, 2) * Power(r, 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                       + Power(a, 2) * Power(r, 4) + Power(r, 6)));
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::l_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {

  auto r = r_KS(a, x, y, z);
  metric_server::spatial_vector lshift{};

  lshift[0] = (2 * M * Power(r, 3) * (a * y + x * r))
              / ((Power(a, 2) + Power(r, 2)) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  lshift[1] = (2 * M * Power(r, 3) * (-(a * x) + y * r))
              / ((Power(a, 2) + Power(r, 2)) * (Power(a, 2) * Power(z, 2) + Power(r, 4)));

  lshift[2] = (2 * M * z * Power(r, 2)) / (Power(a, 2) * Power(z, 2) + Power(r, 4));

  return lshift;
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::u_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {

  auto r = r_KS(a, x, y, z);

  metric_server::spatial_vector ushift{};

  ushift[0] = (2 * M * Power(r, 3) * (a * y + x * r))
              / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                 + Power(a, 2) * Power(r, 4) + Power(r, 6));

  ushift[1] = (2 * M * Power(r, 3) * (-(a * x) + y * r))
              / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                 + Power(a, 2) * Power(r, 4) + Power(r, 6));

  ushift[2] = (2 * M * z * Power(r, 2) * (Power(a, 2) + Power(r, 2)))
              / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r
                 + Power(a, 2) * Power(z, 2) * Power(r, 2)
                 + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r, 3)
                 + Power(a, 2) * Power(r, 4) + Power(r, 6));

  return ushift;
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