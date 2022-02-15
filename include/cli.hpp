#ifndef GRLENSING_CLI_HPP
#define GRLENSING_CLI_HPP

#include "kernel.hpp"
#include "time_integration.hpp"

#include <filesystem>
#include <string>
#include <vector>
#include <yaml-cpp/node/node.h>

namespace grlensing {

/**
 * The program help/usage message that is also used to generate the command line parser.
 */
static constexpr const char *USAGE =
    R"(GRLensing.

    Usage:
      grlensing list-plugins
      grlensing dump-metric <metric-name> [--extension=<ext>]
      grlensing integrate-trajectory <trajectory-file>
      grlensing penrose-breakup <breakup-file>
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
void push_to_paths(const std::filesystem::path &plugin_folder, const YAML::Node &array,
                   const std::string &suffix, std::vector<std::filesystem::path> &path_vector);

/**
 * Construct the path objects of all plugins to be loaded.
 *
 * @param config_file The general yalm configuration file
 * @return A vector containig the path of all plugins to be loaded.
 */
auto make_plugin_paths(const YAML::Node &config_file) -> std::vector<std::filesystem::path>;

/**
 * Loads all plugins in the configured paths.
 *
 * @param kernel The application kernel
 * @param config_file The yaml general configuration file.
 */
void load_plugins(grlensing::kernel &kernel, const YAML::Node &config_file);

/**
 * Implements the list-plugins functionality.
 *
 * The list-plugins mode is usefull for finding out how many plugins and what types of plugins were
 * loaded and are available to be used during execution.
 * TODO: Make this take a const  reference to a kernel
 *
 * @param kernel The application kernel.
 */
void list_plugins(grlensing::kernel &kernel);

/**
 * Implements the integrate-trajectory mode.
 *
 * This mode integrates a single trajectory with parameters specified in a yaml configuration file.
 * This is usefull for doing single trajectory tests to quickly visualize or do physics with a
 * single trajectory.
 * TODO: Make this take a cosnt kernel ref.
 */
void integrate_trajectory(grlensing::kernel &kernel, const grlensing::integrator_config &int_conf,
                          const std::string &traj_file);

/**
 * Implements the dump-metric command.
 *
 * The dump-metric command is usefull for testing spacetime metric implementations. All metric
 * quantities are output as per the confoguration in the general yaml configuration file.
 *
 * @param kernel The application kernel.
 * @param metric_name The name of the spacetime metric to dump.
 * @param config_file The general yaml config file.
 * TODO: Take a const ref to a kernel
 */
void dump_metric(grlensing::kernel &kernel, const std::string &metric_name,
                 const YAML::Node &config_file);

/**
 * Implements the penrose-breakup mode.
 *
 * TODO: Doc.
 * TODO: Make this take a cosnt kernel ref.
 */
void penrose_breakup(grlensing::kernel &kernel, const grlensing::integrator_config &int_conf,
                     const std::string &break_file);

} // namespace grlensing

#endif // GRLENSING_CLI_HPP