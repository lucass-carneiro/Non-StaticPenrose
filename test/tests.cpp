#include <gtest/gtest.h>
#include <mpi.h>

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
  MPI::Init(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  auto test_run_code = RUN_ALL_TESTS();

  MPI::Finalize();
  return test_run_code;
}