#include "cli.hpp"
#include "log.hpp"

#include <docopt/docopt.h>
#include <yaml-cpp/exceptions.h>

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

    } else if (args["penrose-breakup"].asBool()) {
      penrose_breakup(proc_kernel, int_conf, args["<breakup-file>"].asString());
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
  } catch (YAML::Exception &e) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      log<LogEvent::error>("Yaml parser error: {:s}", e.what());
    }
  } catch (std::exception &e) {
    fmt::print("Error: {:s}\n", e.what());
  }

  MPI::Finalize();
  return 0;
}
