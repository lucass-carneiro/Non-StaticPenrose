#include "cli.hpp"

void grlensing::integrate_trajectory(grlensing::kernel &kernel,
                                     const grlensing::integrator_config &int_conf,
                                     const std::string &traj_file) {
  using grlensing::log;
  using grlensing::LogEvent;
  using grlensing::trajectory_config;

  log<LogEvent::message>("Single trajectory mode started. Trajectory  file: {:s}", traj_file);

  // Open trajectory file and load trajectory config
  const auto trajectory_file = YAML::LoadFile(traj_file);
  trajectory_config traj_conf(trajectory_file);

  const auto &output_file = std::filesystem::path(trajectory_file["output_file"].as<std::string>());
  const auto &metadata_file
      = std::filesystem::path(trajectory_file["metadata_file"].as<std::string>());

  // Find a writer capable of outputing the requested data file.
  const auto &writer = kernel.get_storage_server().get_writer(output_file.extension().string());

  // Find a metric with the requested name.
  const auto &metric = kernel.get_metric_server().get_metric(
      trajectory_file["background_metric"].as<std::string>());

  // Load metric parameters
  metric->load_parameters(trajectory_file);

  // Normalize the velocities to ensure that they are physical
  normalize(traj_conf, metric);

  // Open output files for writing.
  writer->open_file(output_file, metadata_file);

  // Integrate trajectory
  integrate(int_conf, writer, metric, traj_conf);

  // Close oppened files
  writer->close_file();
}