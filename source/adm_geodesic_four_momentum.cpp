#include "metric_server.hpp"
#include "time_integration.hpp"

auto grlensing::reconstruct_u_p(const metric_server::metric_ptr &metric, double ti,
                                const metric_server::spatial_vector &Vi,
                                const metric_server::spatial_vector &Xi, double En)
    -> std::array<double, 4> {

  const auto lapse = metric->lapse(ti, Xi[0], Xi[1], Xi[2]);
  const auto u_shift = metric->u_shift(ti, Xi[0], Xi[1], Xi[2]);

  std::array<double, 4> pmu = {En / lapse, Vi[0] - u_shift[0] / lapse, Vi[1] - u_shift[1] / lapse,
                               Vi[2] - u_shift[2] / lapse};

  return pmu;
}