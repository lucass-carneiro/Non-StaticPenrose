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
  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  return Sqrt(1
              - (2 * M * Power(r(x, y, z), 3) * (Power(a, 2) + Power(r(x, y, z), 2)))
                    / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
                       + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
                       + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))); // NOLINT
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::l_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {
  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  return metric_server::spatial_vector{(2 * M * Power(r(x, y, z), 3) * (a * y + x * r(x, y, z)))
                                           / ((Power(a, 2) + Power(r(x, y, z), 2))
                                              * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))),
                                       (2 * M * Power(r(x, y, z), 3) * (-(a * x) + y * r(x, y, z)))
                                           / ((Power(a, 2) + Power(r(x, y, z), 2))
                                              * (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))),
                                       (2 * M * z * Power(r(x, y, z), 2))
                                           / (Power(a, 2) * Power(z, 2) + Power(r(x, y, z), 4))};
}

GRLENSING_KERRSCHILD_KERR_METRIC_API auto KerrSchild_Kerr::u_shift(double, double x, double y,
                                                                   double z)
    -> metric_server::spatial_vector {
  auto r = [&](double x, double y, double z) { return r_KS(a, x, y, z); };
  return metric_server::spatial_vector{
      (2 * M * Power(r(x, y, z), 3) * (a * y + x * r(x, y, z)))
          / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
             + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
             + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
             + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6)),
      (2 * M * Power(r(x, y, z), 3) * (-(a * x) + y * r(x, y, z)))
          / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
             + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
             + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
             + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6)),
      (2 * M * z * Power(r(x, y, z), 2) * (Power(a, 2) + Power(r(x, y, z), 2)))
          / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * r(x, y, z)
             + Power(a, 2) * Power(z, 2) * Power(r(x, y, z), 2)
             + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(r(x, y, z), 3)
             + Power(a, 2) * Power(r(x, y, z), 4) + Power(r(x, y, z), 6))};
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