
#include "isotropic_Schwarzschild.hpp"

using grlensing::isotropic_Schwarzschild;
using grlensing::kernel;
using grlensing::metric_server;

template <typename T> constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T(1) : x * Power(x, n - 1);
}

template <typename T> constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T(1) / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

constexpr auto Sqrt(auto x) { return std::sqrt(x); }

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API void
isotropic_Schwarzschild::load_parameters(const YAML::Node &node) {
  M = node["Isotropic_Schwarzschild_Settings"]["M"].as<double>();
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto isotropic_Schwarzschild::lapse(double, double x,
                                                                                 double y, double z)
    -> double {
  return 1 - (2 * M) / (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2)));
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto isotropic_Schwarzschild::l_shift(double, double,
                                                                                   double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{0.0, 0.0, 0.0};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto isotropic_Schwarzschild::u_shift(double, double,
                                                                                   double, double)
    -> metric_server::spatial_vector {
  return metric_server::spatial_vector{0.0, 0.0, 0.0};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::ll_smetric(double, double x, double y, double z)
    -> metric_server::spatial_matrix {
  const auto nonzero_term = Power(1 + M / (2. * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))), 4);
  return metric_server::spatial_matrix{
      {{nonzero_term, 0.0, 0.0}, {0.0, nonzero_term, 0.0}, {0.0, 0.0, nonzero_term}}};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::uu_smetric(double, double x, double y, double z)
    -> metric_server::spatial_matrix {
  const auto nonzero_term
      = 1.0 / (Power(1 + M / (2. * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))), 4));
  return metric_server::spatial_matrix{
      {{nonzero_term, 0.0, 0.0}, {0.0, nonzero_term, 0.0}, {0.0, 0.0, nonzero_term}}};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::ll_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{{{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::ul_extrinsic(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{{{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::spatial_christoffel(double, double x, double y, double z)
    -> metric_server::chirstofell_t {
  metric_server::chirstofell_t Gamma;

  Gamma[0][0][0] = (-2 * M * x)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[0][0][1] = (-2 * M * y)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[0][0][2] = (-2 * M * z)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[0][1][1] = (2 * M * x)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[0][1][2] = 0.0;
  Gamma[0][2][2] = (2 * M * x)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));

  Gamma[1][0][0] = (2 * M * y)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[1][0][1] = (-2 * M * x)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[1][0][2] = 0.0;
  Gamma[1][1][1] = (-2 * M * y)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[1][1][2] = (-2 * M * z)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[1][2][2] = (2 * M * y)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));

  Gamma[2][0][0] = (2 * M * z)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[2][0][1] = 0.0;
  Gamma[2][0][2] = (-2 * M * x)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[2][1][1] = (2 * M * z)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[2][1][2] = (-2 * M * y)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));
  Gamma[2][2][2] = (-2 * M * z)
                   / ((Power(x, 2) + Power(y, 2) + Power(z, 2))
                      * (M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))));

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

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::grad_lapse(double, double x, double y, double z)
    -> metric_server::spatial_vector {

  metric_server::spatial_vector dalpha;

  dalpha[0] = (4 * M * x)
              / (Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))
                 * Power(M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2)), 2));
  dalpha[1] = (4 * M * y)
              / (Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))
                 * Power(M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2)), 2));
  dalpha[2] = (4 * M * z)
              / (Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))
                 * Power(M + 2 * Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2)), 2));

  return dalpha;
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto
isotropic_Schwarzschild::grad_ushift(double, double, double, double)
    -> metric_server::spatial_matrix {
  return metric_server::spatial_matrix{{{0.0, 0.0, 0.0}, {0.0, 0.0, 0}, {0.0, 0.0, 0.0}}};
}

GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto isotropic_Schwarzschild::name()
    -> std::string_view {
  return "Isotropic Schwarzschild";
}

extern "C" GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API auto get_engine_version() -> unsigned {
  return unsigned(1);
}

extern "C" GRLENSING_ISOTROPIC_SCHWARZSCHILD_METRIC_API void register_plugin(kernel &kernel) {
  kernel.get_metric_server().add_metric(
      std::unique_ptr<metric_server::adm_metric>(new isotropic_Schwarzschild()));
}
