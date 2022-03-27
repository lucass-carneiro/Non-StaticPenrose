#ifndef GRLENSING_TESTS_TEST_CONSTANTS_HPP
#define GRLENSING_TESTS_TEST_CONSTANTS_HPP

#include <cmath>

namespace grlensing_tests {

/**
 * TODO: doc
 */
constexpr const double double_comp_tol = 1.0e-15;

/**
 * TODO: doc
 */
constexpr const int random_seed = 100;

/**
 * TODO: doc
 */
constexpr const double coord_range[2] = {-10.0, 10.0};

/**
 * TODO: doc
 */
constexpr const double par_range[2] = {0.1, 1.0};

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

} // namespace grlensing_tests

#endif // GRLENSING_TESTS_TEST_CONSTANTS_HPP