#include "metric_server.hpp"
#include "time_integration.hpp"

auto grlensing::compute_global_energy(const metric_server::metric_ptr &metric, double ti,
                                      const std::array<double, 3> &Vi,
                                      const std::array<double, 3> &Xi, double El) -> double {

  const auto lapse = metric->lapse(ti, Xi[0], Xi[1], Xi[2]);
  const auto u_shift = metric->u_shift(ti, Xi[0], Xi[1], Xi[2]);
  const auto ll_smetric = metric->ll_smetric(ti, Xi[0], Xi[1], Xi[2]);

  // $\gamma_{ij} \beta^i V^j
  auto contraction
      = ll_smetric[0][0] * u_shift[0] * Vi[0] + ll_smetric[1][1] * u_shift[1] * Vi[1]
        + ll_smetric[2][2] * u_shift[2] * Vi[2] + 2 * ll_smetric[0][1] * u_shift[0] * Vi[1]
        + 2 * ll_smetric[0][2] * u_shift[0] * Vi[2] + 2 * ll_smetric[1][2] * u_shift[1] * Vi[2];

  return (lapse - contraction) * El;
}