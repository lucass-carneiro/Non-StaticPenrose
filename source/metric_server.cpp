#include "metric_server.hpp"

#include "log.hpp"

GRLENSING_API void grlensing::metric_server::print_names() const {
  for (const auto &metric : metrics) {
    log<LogEvent::message>("Spacetime metric: {:s}", metric->name());
  }
}

GRLENSING_API auto grlensing::metric_server::get_metric(std::string_view name) const noexcept(false)
    -> const std::unique_ptr<grlensing::metric_server::adm_metric> & {
  for (const auto &metric : metrics) {
    if (metric->name() == name)
      return metric;
  }
  log<LogEvent::error>("No plugin available that provides the metric with name {:s}", name);
  throw std::runtime_error("IO error");
}
auto grlensing::shift_module(const metric_server::metric_ptr &metric, double t, double x, double y,
                             double z) -> double {
  const auto l_shift = metric->l_shift(t, x, y, z);
  const auto u_shift = metric->u_shift(t, x, y, z);

  return l_shift[0] * u_shift[0] + l_shift[1] * u_shift[1] + l_shift[2] * u_shift[2];
}

auto grlensing::ll_g_00(const metric_server::metric_ptr &metric, double t, double x, double y,
                        double z) -> double {
  const auto shift_mod = shift_module(metric, t, x, y, z);
  const auto lapse = metric->lapse(t, x, y, z);

  return shift_mod - lapse * lapse;
}

auto grlensing::compute_global_energy(const metric_server::metric_ptr &metric, double ti,
                                      const metric_server::spatial_vector &Vi,
                                      const metric_server::spatial_vector &Xi, double El)
    -> double {

  const auto lapse = metric->lapse(ti, Xi[0], Xi[1], Xi[2]);
  const auto u_shift = metric->u_shift(ti, Xi[0], Xi[1], Xi[2]);
  const auto ll_smetric = metric->ll_smetric(ti, Xi[0], Xi[1], Xi[2]);

  // $\gamma_{ij} \beta^i V^j
  auto contraction
      = ll_smetric[0][0] * u_shift[0] * Vi[0] + ll_smetric[1][1] * u_shift[1] * Vi[1]
        + ll_smetric[2][2] * u_shift[2] * Vi[2] + 2 * ll_smetric[0][1] * u_shift[0] * Vi[1]
        + 2 * ll_smetric[0][2] * u_shift[0] * Vi[2] + 2 * ll_smetric[1][2] * u_shift[1] * Vi[2];

  return (lapse - contraction) * El; // TODO: Change this for the non unit mass
}

auto grlensing::reconstruct_u_p(const metric_server::metric_ptr &metric, double ti,
                                const metric_server::spatial_vector &Vi,
                                const metric_server::spatial_vector &Xi, double En)
    -> std::array<double, 4> {

  const auto lapse = metric->lapse(ti, Xi[0], Xi[1], Xi[2]);
  const auto u_shift = metric->u_shift(ti, Xi[0], Xi[1], Xi[2]);

  std::array<double, 4> pmu
      = {En / lapse, En * (Vi[0] - u_shift[0]) / lapse, En * (Vi[1] - u_shift[1]) / lapse,
         En * (Vi[2] - u_shift[2]) / lapse};

  return pmu;
}

auto grlensing::decompose_u_p(const metric_server::metric_ptr &metric,
                              const std::array<double, 4> &u_p, double ti,
                              const metric_server::spatial_vector &Xi)
    -> std::tuple<metric_server::spatial_vector, double> {

  const auto u_shift = metric->u_shift(ti, Xi[0], Xi[1], Xi[2]);
  const auto lapse = metric->lapse(ti, Xi[0], Xi[1], Xi[2]);

  double En = lapse * u_p[0];

  const double denominator = (lapse * u_p[0]);
  metric_server::spatial_vector Vi
      = {(u_p[1] + u_p[0] * u_shift[0]) / denominator, (u_p[2] + u_p[0] * u_shift[1]) / denominator,
         (u_p[3] + u_p[0] * u_shift[2]) / denominator};

  return std::make_tuple(Vi, En);
}