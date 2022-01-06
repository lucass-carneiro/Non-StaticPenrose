#include "cli.hpp"
#include "mpimatrix.hpp"

#include <string>

void grlensing::dump_metric(grlensing::kernel &kernel, const std::string &metric_name,
                            const YAML::Node &config_file) {
  using grlensing::log;
  using grlensing::LogEvent;

  const auto extension = config_file["dump_metric_settings"]["extension"].as<std::string>();
  auto rf = config_file["dump_metric_settings"]["radius"].as<double>();
  auto points = config_file["dump_metric_settings"]["points"].as<unsigned>();

  const auto &writer = kernel.get_storage_server().get_writer(extension);
  const auto &metric = kernel.get_metric_server().get_metric(metric_name);
  const std::filesystem::path lapse_file(metric_name + std::string("_lapse_")
                                         + std::to_string(MPI::COMM_WORLD.Get_rank()) + extension);
  writer->open_file(lapse_file);

  if (points <= 1) {
    log<LogEvent::error>("The number of points in a metric dump must be at least 2");
    throw std::runtime_error("Invalid parameter error");
  }

  auto r0 = -rf;
  const auto dr = 2 * rf / points;

  mpimatrix<double> xy_coords(points + 1, points + 1);
  const auto owned_linear_range = xy_coords.get_owned_linear_range();
  const auto start_ij = xy_coords.global_linear_to_matrix(owned_linear_range.first);
  const auto end_ij = xy_coords.global_linear_to_matrix(owned_linear_range.second);

  log<LogEvent::message>("Dumping metric {:s}", metric_name);

  for (auto i = start_ij.first; i <= end_ij.first; i++) {
    for (auto j = start_ij.second; j <= end_ij.second; j++) {
      const auto x = r0 + i * dr;
      const auto y = r0 + j * dr;

      writer->push_real(x);
      writer->push_real(y);

      writer->push_final_real(metric->lapse(0.0, x, y, 0.0));
    }
  }

  writer->close_file();
}