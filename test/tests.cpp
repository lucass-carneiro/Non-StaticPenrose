#include "kerr_schild_kerr_tests.hpp"

#include <gtest/gtest.h>
#include <mpi.h>
#include <yaml-cpp/exceptions.h>

auto main(int argc, char **argv) -> int {
  MPI::Init(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  int test_run_code = RUN_ALL_TESTS();

  MPI::Finalize();
  return test_run_code;
}