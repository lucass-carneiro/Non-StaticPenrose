#include "cli.hpp"
#include "log.hpp"

#include <gtest/gtest.h>
#include <mpi.h>
#include <yaml-cpp/exceptions.h>

constexpr auto Factorial(unsigned int number) -> unsigned int {
  return (number == 0 || number == 1) ? 1 : number * Factorial(number - 1);
}

TEST(FactorialTest, HandlesZeroInput) { EXPECT_EQ(Factorial(0), 1); }

// Tests factorial of positive numbers.
TEST(FactorialTest, HandlesPositiveInput) {
  EXPECT_EQ(Factorial(1), 1);
  EXPECT_EQ(Factorial(2), 2);
  EXPECT_EQ(Factorial(3), 6);
  EXPECT_EQ(Factorial(8), 40320);
}

auto main(int argc, char **argv) -> int {
  using namespace grlensing;

  int test_run_code = -1;

  try {
    MPI::Init(argc, argv);

    // Load global configuration file
    const auto config_file = YAML::LoadFile("../../configs/grlensing_testing_config.yaml");

    // Initilize kernel
    kernel proc_kernel;

    // Load plugins
    load_plugins(proc_kernel, config_file);

    // Integrator base configuration
    const integrator_config int_conf(config_file);

    // Perform tests
    ::testing::InitGoogleTest(&argc, argv);
    test_run_code = RUN_ALL_TESTS();

  } catch (YAML::Exception &e) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      log<LogEvent::error>("Yaml parser error: {:s}", e.what());
    }
  } catch (std::exception &e) {
    fmt::print("Error: {:s}\n", e.what());
  }

  MPI::Finalize();
  return test_run_code;
}