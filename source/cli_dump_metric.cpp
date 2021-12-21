#include "cli.hpp"

void grlensing::dump_metric(grlensing::kernel &kernel, const std::string &metric_name,
                            const YAML::Node &config_file) {
  using grlensing::log;
  using grlensing::LogEvent;

  const auto extension = config_file["dump_metric_settings"]["extension"].as<std::string>();
  auto rf = config_file["dump_metric_settings"]["radius"].as<double>();
  auto points = config_file["dump_metric_settings"]["points"].as<unsigned>();

  const auto &writer = kernel.get_storage_server().get_writer(extension);
  const auto &metric = kernel.get_metric_server().get_metric(metric_name);
  const std::filesystem::path lapse_file(metric_name + std::string("_lapse") + extension);
  writer->open_file(lapse_file);

  if (points <= 1) {
    log<LogEvent::error>("The number of points in a metric dump must be at least 2");
    throw std::runtime_error("Invalid parameter error");
  } else {
    points--;
  }

  auto r0 = -rf;
  const auto dr = (rf - r0) / points;

  log<LogEvent::message>("Dumping metric {:s}", metric_name);

  for (unsigned i = 0; i < points + 1; i++) {
    for (unsigned j = 0; j < points + 1; j++) {
      for (unsigned k = 0; k < points + 1; k++) {
        const auto x = r0 + i * dr;
        const auto y = r0 + j * dr;
        const auto z = r0 + k * dr;

        writer->push_real(x);
        writer->push_real(y);
        writer->push_real(z);

        writer->push_final_real(metric->lapse(0, x, y, z));
      }
    }
  }

  writer->close_file();
}