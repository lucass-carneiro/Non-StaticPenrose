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