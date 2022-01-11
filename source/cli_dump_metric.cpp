#include "cli.hpp"
#include "mpi_index_map_3D.hpp"

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

  mpi_index_map_3D coords_map(points + 1, points + 1, points + 1);
  const auto glor = coords_map.get_glor();
  const auto start_ijk = coords_map.global_linear_to_matrix(glor.first);
  const auto end_ijk = coords_map.global_linear_to_matrix(glor.second);

  log<LogEvent::message>("Dumping metric {:s}", metric_name);

  for (auto i = std::get<0>(start_ijk); i <= std::get<0>(end_ijk); i++) {
    for (auto j = std::get<1>(start_ijk); j <= std::get<1>(end_ijk); j++) {
      for (auto k = std::get<2>(start_ijk); k <= std::get<2>(end_ijk); k++) {
        const auto x = r0 + i * dr;
        const auto y = r0 + j * dr;
        const auto z = r0 + k * dr;

        writer->push_real(x);
        writer->push_real(y);
        writer->push_real(z);

        writer->push_final_real(metric->lapse(0.0, x, y, z));
      }
    }
  }

  writer->close_file();
}