#include "cli.hpp"
#include "log.hpp"
#include "metric_server.hpp"
#include "time_integration.hpp"

void grlensing::penrose_breakup(grlensing::kernel &kernel,
                                const grlensing::integrator_config &int_conf,
                                const std::string &break_file) {

  using grlensing::log;
  using grlensing::LogEvent;
  using grlensing::trajectory_config;

  log<LogEvent::message>("Penrose particle breakup mode started. Breakup  file: {:s}", break_file);

  // Open breakup file and load trajectory configs
  const auto breakup_file = YAML::LoadFile(break_file);

  // Find a writer capable of outputing the requested data file.
  const auto &writer
      = kernel.get_storage_server().get_writer(breakup_file["output_format"].as<std::string>());

  // Find a metric with the requested name.
  const auto &metric
      = kernel.get_metric_server().get_metric(breakup_file["background_metric"].as<std::string>());

  // Load metric parameters
  metric->load_parameters(breakup_file);

  // Trajectory 1 file name and path
  const auto traj_1_file_name
      = std::string("trajectory_1") + breakup_file["output_format"].as<std::string>();
  const auto traj_1_metadata_file_name
      = std::string("trajectory_1_metadata") + breakup_file["output_format"].as<std::string>();
  const auto traj_1_file = std::filesystem::path(traj_1_file_name);
  const auto traj_1_metadata_file = std::filesystem::path(traj_1_metadata_file_name);

  // Trajectory 2 file name and path
  const auto traj_2_file_name
      = std::string("trajectory_2") + breakup_file["output_format"].as<std::string>();
  const auto traj_2_metadata_file_name
      = std::string("trajectory_2_metadata") + breakup_file["output_format"].as<std::string>();
  const auto traj_2_file = std::filesystem::path(traj_2_file_name);
  const auto traj_2_metadata_file = std::filesystem::path(traj_2_metadata_file_name);

  // Trajectory 3 file name and path
  const auto traj_3_file_name
      = std::string("trajectory_3") + breakup_file["output_format"].as<std::string>();
  const auto traj_3_metadata_file_name
      = std::string("trajectory_3_metadata") + breakup_file["output_format"].as<std::string>();
  const auto traj_3_file = std::filesystem::path(traj_3_file_name);
  const auto traj_3_metadata_file = std::filesystem::path(traj_3_metadata_file_name);

  /* The spatial vector containing the breakup point. This is the initial position common to all
   * trajectories
   */
  metric_server::spatial_vector breakup_point
      = {breakup_file["breakup_point"]["X1"].as<realtype>(),
         breakup_file["breakup_point"]["X2"].as<realtype>(),
         breakup_file["breakup_point"]["X3"].as<realtype>()};

  // First trajectory configuration
  metric_server::spatial_vector V_traj1
      = {breakup_file["trajectory_1"]["initial_V1"].as<realtype>(),
         breakup_file["trajectory_1"]["initial_V2"].as<realtype>(),
         breakup_file["trajectory_1"]["initial_V3"].as<realtype>()};

  trajectory_config traj_1_conf;
  traj_1_conf.initial_time = 0.0;
  traj_1_conf.final_time = breakup_file["trajectory_1"]["final_time"].as<realtype>();
  traj_1_conf.initial_V1 = V_traj1[0];
  traj_1_conf.initial_V2 = V_traj1[1];
  traj_1_conf.initial_V3 = V_traj1[2];
  traj_1_conf.initial_X1 = breakup_point[0];
  traj_1_conf.initial_X2 = breakup_point[1];
  traj_1_conf.initial_X3 = breakup_point[2];
  traj_1_conf.initial_EN = breakup_file["trajectory_1"]["initial_EN"].as<realtype>();
  traj_1_conf.background_radius = breakup_file["background_radius"].as<realtype>();
  traj_1_conf.energy_threshold = breakup_file["energy_threshold"].as<realtype>();
  traj_1_conf.output_times = breakup_file["trajectory_1"]["output_times"].as<int>();
  traj_1_conf.particle_type = breakup_file["particle_type"].as<int>();

  // Second trajectory configuration
  metric_server::spatial_vector V_traj2
      = {breakup_file["trajectory_2"]["initial_V1"].as<realtype>(),
         breakup_file["trajectory_2"]["initial_V2"].as<realtype>(),
         breakup_file["trajectory_2"]["initial_V3"].as<realtype>()};

  trajectory_config traj_2_conf;
  traj_2_conf.initial_time = 0.0;
  traj_2_conf.final_time = breakup_file["trajectory_2"]["final_time"].as<realtype>();
  traj_2_conf.initial_V1 = V_traj2[0];
  traj_2_conf.initial_V2 = V_traj2[1];
  traj_2_conf.initial_V3 = V_traj2[2];
  traj_2_conf.initial_X1 = breakup_point[0];
  traj_2_conf.initial_X2 = breakup_point[1];
  traj_2_conf.initial_X3 = breakup_point[2];
  traj_2_conf.initial_EN = breakup_file["trajectory_2"]["initial_EN"].as<realtype>();
  traj_2_conf.background_radius = breakup_file["background_radius"].as<realtype>();
  traj_2_conf.energy_threshold = breakup_file["energy_threshold"].as<realtype>();
  traj_2_conf.output_times = breakup_file["trajectory_2"]["output_times"].as<int>();
  traj_2_conf.particle_type = breakup_file["particle_type"].as<int>();

  // Apply 4-momentum conservation
  auto p_traj1 = reconstruct_u_p(metric, 0.0, V_traj1, breakup_point, traj_1_conf.initial_EN);
  auto p_traj2 = reconstruct_u_p(metric, 0.0, V_traj2, breakup_point, traj_2_conf.initial_EN);
  std::array<double, 4> p_traj3 = {p_traj1[0] - p_traj2[0], p_traj1[1] - p_traj2[1],
                                   p_traj1[2] - p_traj2[2], p_traj1[3] - p_traj2[3]};

  // Using p_traj3, we solve for EN, V1, V3 and V3 for the final trajectory
  auto conserved_data = decompose_u_p(metric, p_traj3, 0.0, breakup_point);

  trajectory_config traj_3_conf;
  traj_3_conf.initial_time = 0.0;
  traj_3_conf.final_time = breakup_file["trajectory_3"]["final_time"].as<realtype>();
  traj_3_conf.initial_V1 = std::get<0>(conserved_data)[0];
  traj_3_conf.initial_V2 = std::get<0>(conserved_data)[1];
  traj_3_conf.initial_V3 = std::get<0>(conserved_data)[2];
  traj_3_conf.initial_X1 = breakup_point[0];
  traj_3_conf.initial_X2 = breakup_point[1];
  traj_3_conf.initial_X3 = breakup_point[2];
  traj_3_conf.initial_EN = std::get<1>(conserved_data);
  traj_3_conf.background_radius = breakup_file["background_radius"].as<realtype>();
  traj_3_conf.energy_threshold = breakup_file["energy_threshold"].as<realtype>();
  traj_3_conf.output_times = breakup_file["trajectory_3"]["output_times"].as<int>();
  traj_3_conf.particle_type = breakup_file["particle_type"].as<int>();

  metric_server::spatial_vector V_traj3
      = {traj_3_conf.initial_V1, traj_3_conf.initial_V2, traj_3_conf.initial_V3};

  // Normalize velocities
  normalize(traj_1_conf, metric);
  normalize(traj_2_conf, metric);
  normalize(traj_3_conf, metric);

  // Integrate trajectory 1
  writer->open_file(traj_1_file, traj_1_metadata_file);

  integrate(int_conf, writer, metric, traj_1_conf);

  const auto mass1 = compute_mass(metric, 0, V_traj1, breakup_point, traj_1_conf.initial_EN);
  writer->push_metadata("   mass", mass1);

  writer->close_file();

  // Integrate trajectory 2
  writer->open_file(traj_2_file, traj_2_metadata_file);

  integrate(int_conf, writer, metric, traj_2_conf);

  const auto mass2 = compute_mass(metric, 0, V_traj2, breakup_point, traj_2_conf.initial_EN);
  writer->push_metadata("   mass", mass2);

  writer->close_file();

  // Integrate trajectory 3
  writer->open_file(traj_3_file, traj_3_metadata_file);

  integrate(int_conf, writer, metric, traj_3_conf);

  const auto mass3 = compute_mass(metric, 0, V_traj3, breakup_point, traj_3_conf.initial_EN);
  writer->push_metadata("   mass", mass3);

  writer->close_file();
}