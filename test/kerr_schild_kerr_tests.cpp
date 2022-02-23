#include "kerr_schild_kerr_tests.hpp"

#include "cli.hpp"
#include "log.hpp"
#include "test_constants.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <mpi.h>
#include <random>

template <typename T> inline constexpr auto Power(T x, unsigned n) -> double {
  return (n == 0) ? T{1} : x * Power(x, n - 1);
}

template <typename T> inline constexpr auto Power(T x, int n) -> double {
  return (n < 0) ? T{1} / Power(x, unsigned(-n)) : Power(x, unsigned(n));
}

inline constexpr auto Sqrt(auto x) {
  using std::sqrt;
  return sqrt(x);
}

auto grlensing_tests::r(double M, double a, double x, double y, double z) -> double {
  using std::sqrt;
  double part1 = -a * a + x * x + y * y + z * z;
  return sqrt(0.5 * (part1 + sqrt(4 * a * a * z * z + part1 * part1)));
}

auto grlensing_tests::alpha(double M, double a, double x, double y, double z) -> double {
  double rVal = r(M, a, x, y, z);
  return Sqrt(1
              - (2 * M * Power(rVal, 3) * (Power(a, 2) + Power(rVal, 2)))
                    / (Power(a, 4) * Power(z, 2) + 2 * Power(a, 2) * M * Power(z, 2) * rVal
                       + Power(a, 2) * Power(z, 2) * Power(rVal, 2)
                       + 2 * M * (Power(x, 2) + Power(y, 2) + Power(z, 2)) * Power(rVal, 3)
                       + Power(a, 2) * Power(rVal, 4) + Power(rVal, 6)));
}

TEST(Kerr_Schild_Kerr_Tests, correct_lapse) {
  using namespace grlensing_tests;
  using namespace grlensing;
  using std::mt19937_64;
  using std::uniform_real_distribution;

  try {
    // Load global configuration file
    const auto config_file = YAML::LoadFile("../../configs/grlensing_testing_config.yaml");

    // Initilize kernel
    kernel proc_kernel;

    // Load plugins
    load_plugins(proc_kernel, config_file);

    // Random parameter generation
    mt19937_64 engine(random_seed);
    uniform_real_distribution<double> mass_distrib(par_range[0], par_range[1]);
    uniform_real_distribution<double> coord_distrib(coord_range[0], coord_range[1]);

    const double M = 1;
    const double a = M / 2;

    const double x = coord_distrib(engine);
    const double y = coord_distrib(engine);
    const double z = coord_distrib(engine);

    const auto &metric = proc_kernel.get_metric_server().get_metric("Kerr-Schild Kerr");

    // Tests
    EXPECT_NEAR(alpha(M, a, x, y, z), metric->lapse(0.0, x, y, z), double_comp_tol);

  } catch (YAML::Exception &e) {
    if (MPI::COMM_WORLD.Get_rank() == 0) {
      log<LogEvent::error>("Yaml parser error: {:s}", e.what());
    }
  } catch (std::exception &e) {
    fmt::print("Error: {:s}\n", e.what());
  }
}