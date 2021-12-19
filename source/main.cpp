#include "kernel.hpp"
#include "log.hpp"
#include "options.hpp"
#include "plugin.hpp"
#include "time_integration.hpp"

#include <cstddef>
#include <docopt.h>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sundials/sundials_types.h>
#include <vector>
#include <yaml-cpp/node/node.h>

namespace fs = std::filesystem;

static constexpr const char *USAGE =
    R"(GRLensing.

    Usage:
      grlensing list-plugins
      grlensing dump-metric <metric-name> [--extension=<ext>]
      grlensing integrate-trajectory <trajectory-file>
      grlensing (-h | --help)
      grlensing --version

    Options:
      -h --help          Show this screen.
      --version          Show version.
      --extension=<ext>  The extension of the output file. [default: .ascii]
)";

/**
 * Adds a path object to an array contaning all the paths of the plugins to be loaded.
 *
 * @param plugin_folder A path to the base folder where plugins are stored.
 * @param array A node to a Yaml file containing an array of plugin names to load.
 * @param suffix The suffix that must be appended to the name of the library files.
 * @param path_vector A vector with all the paths of plugins to be loaded.
 */
static void push_to_paths(const fs::path &plugin_folder, const YAML::Node &array,
                          const std::string &suffix, std::vector<fs::path> &path_vector) {
  for (const auto &plugin_name : array) {
    const auto base_name = plugin_name.as<std::string>();

    // TODO: Find a better way to do this
    const auto lib_name = base_name + std::string("/libgrlensing_") + base_name + suffix;

    const auto lib_path = fs::path(lib_name);
    path_vector.push_back(plugin_folder / lib_path);
  }
}

static auto make_plugin_paths(const YAML::Node &config_file) -> std::vector<fs::path> {
  // Base plugin folders
  const auto metric_plugin_folder = fs::path(config_file["metric_plugin_folder"].as<std::string>());
  const auto archive_reader_plugin_folder
      = fs::path(config_file["archive_reader_plugin_folder"].as<std::string>());
  const auto archive_writer_plugin_folder
      = fs::path(config_file["archive_writer_plugin_folder"].as<std::string>());

  // Find all plugins in the configured paths
  std::vector<fs::path> plugin_paths;
  plugin_paths.reserve(config_file["load_metrics"].size() + config_file["load_readers"].size());

  push_to_paths(metric_plugin_folder, config_file["load_metrics"], "_metric.so", plugin_paths);
  push_to_paths(archive_reader_plugin_folder, config_file["load_readers"], "_archive.so",
                plugin_paths);
  push_to_paths(archive_writer_plugin_folder, config_file["load_writers"], "_writer.so",
                plugin_paths);

  return plugin_paths;
}

/**
 * TODO: Document.
 */
static void load_plugins(grlensing::kernel &kernel, const YAML::Node &config_file) {
  const auto plugin_paths = make_plugin_paths(config_file);

  for (const auto &plugin_path : plugin_paths) {
    kernel.load_plugin(plugin_path);
  }
}

/**
 * TODO: Document. make this take a const  reference to a kernel
 */
static void list_plugins(grlensing::kernel &kernel) {
  using grlensing::log;
  using grlensing::LogEvent;

  log<LogEvent::message>("Loaded plugins: {:d}", kernel.get_plugin_count());
  kernel.get_storage_server().print_names();
  kernel.get_metric_server().print_names();
}

/**
 * TODO: Doc. Make this take a cosnt kernel ref.
 */
static void integrate_trajectory(grlensing::kernel &kernel,
                                 const grlensing::integrator_config &int_conf,
                                 const std::string &traj_file) {
  using grlensing::log;
  using grlensing::LogEvent;
  using grlensing::trajectory_config;

  log<LogEvent::message>("Single trajectory mode started. Trajectory  file: {:s}", traj_file);

  // Open trajectory file and load trajectory config
  const auto trajectory_file = YAML::LoadFile(traj_file);
  trajectory_config traj_conf(trajectory_file);

  const auto &output_file = fs::path(trajectory_file["output_file"].as<std::string>());
  const auto &metadata_file = fs::path(trajectory_file["metadata_file"].as<std::string>());

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

/**
 * TODO: doc. Take a const ref to a kernel
 * TODO: Finish. Dump metric
 */
void dump_metric(grlensing::kernel &kernel, const std::string &metric_name,
                 const YAML::Node &config_file) {
  using grlensing::log;
  using grlensing::LogEvent;

  const auto extension = config_file["dump_metric_settings"]["extension"].as<std::string>();
  auto rf = config_file["dump_metric_settings"]["radius"].as<double>();
  auto points = config_file["dump_metric_settings"]["points"].as<unsigned>();

  const auto &writer = kernel.get_storage_server().get_writer(extension);
  const auto &metric = kernel.get_metric_server().get_metric(metric_name);
  const fs::path lapse_file(metric_name + std::string("_lapse") + extension);
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

auto main(int argc, char **argv) -> int {
  using namespace grlensing;

  try {
    MPI::Init(argc, argv);

    // Parse command line arguments
    auto args = docopt::docopt_parse(USAGE, {argv + 1, argv + argc}); // NOLINT

    // Load global configuration file
    const auto config_file = YAML::LoadFile("grlensing_config.yaml");

    // Initilize kernel
    kernel proc_kernel;

    // Load plugins
    load_plugins(proc_kernel, config_file);

    // Integrator base configuration
    const integrator_config int_conf(config_file);

    // Take actions for program options
    if (args["list-plugins"].asBool()) {
      list_plugins(proc_kernel);

    } else if (args["dump-metric"].asBool()) {
      dump_metric(proc_kernel, args["<metric-name>"].asString(), config_file);

    } else if (args["integrate-trajectory"].asBool()) {
      integrate_trajectory(proc_kernel, int_conf, args["<trajectory-file>"].asString());
    }
  } catch (docopt::DocoptLanguageError &) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      fmt::print("Internal problem: a syntax error ocurred in the USAGE string. Please contact a "
                 "developper\n");
    }
  } catch (docopt::DocoptExitHelp &) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      fmt::print("{:s}\n", USAGE);
    }
  } catch (docopt::DocoptExitVersion &) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      fmt::print("GRlensing version {:d}.{:d}.{:d}\n", GRLENSING_VERSION_MAJOR,
                 GRLENSING_VERSION_MINOR, GRLENSING_VERSION_PATCH);
    }
  } catch (docopt::DocoptArgumentError &) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      fmt::print(
          "Unrecognized arguments passed. Rerun with the --help option for usage instructions.\n");
    }
  } catch (std::exception &e) {
    fmt::print("Error: {:s}\n", e.what());
  }

  MPI::Finalize();
  return 0;
}
