#include "metric_server.hpp"
#include "time_integration.hpp"

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