#ifndef GRLENSING_TESTS_TEST_CONSTANTS_HPP
#define GRLENSING_TESTS_TEST_CONSTANTS_HPP

#include <cmath>
#include <concepts>
#include <functional>

namespace grlensing_tests {

/**
 * TODO: doc
 */
constexpr const double double_comp_tol = 1.0e-15;

/**
 * The comaprison tolerance when comparing data obtained by finite diferences
 */
constexpr const double fd_comp_tol = 1.0e-8;

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

/**
 * What direction (or function parameter) to use when performin finite differences
 */
enum class fd_direction { t, x, y, z };

/**
 * Computes the fourth order accurate finite difference derivative of a univariate function.
 *
 * @param function The function to derivate
 * @param point The point at which the derivative is to be computed.
 * @param fd_delta The "grid spacing" to use for the finite difference operator.
 * @return The value of the first derivative at the specified point
 */
template <fd_direction dir, std::floating_point T>
inline auto fd_4(std::function<T(T, T, T, T)> function, T t, T x, T y, T z, T fd_delta = T{1.0e-8})
    -> T {
  T point_p_1d = T{0};
  T point_p_2d = T{0};
  T point_m_1d = T{0};
  T point_m_2d = T{0};

  auto f_p_1d = T{0};
  auto f_p_2d = T{0};
  auto f_m_1d = T{0};
  auto f_m_2d = T{0};

  if constexpr (dir == fd_direction::t) {
    point_p_1d = t + fd_delta;
    point_p_2d = t + T{2} * fd_delta;
    point_m_1d = t - fd_delta;
    point_m_2d = t - T{2} * fd_delta;

    f_p_1d = function(point_p_1d, x, y, z);
    f_p_2d = function(point_p_2d, x, y, z);
    f_m_1d = function(point_m_1d, x, y, z);
    f_m_2d = function(point_m_2d, x, y, z);
  } else if constexpr (dir == fd_direction::x) {
    point_p_1d = x + fd_delta;
    point_p_2d = x + T{2} * fd_delta;
    point_m_1d = x - fd_delta;
    point_m_2d = x - T{2} * fd_delta;

    f_p_1d = function(t, point_p_1d, y, z);
    f_p_2d = function(t, point_p_2d, y, z);
    f_m_1d = function(t, point_m_1d, y, z);
    f_m_2d = function(t, point_m_2d, y, z);
  } else if constexpr (dir == fd_direction::y) {
    point_p_1d = y + fd_delta;
    point_p_2d = y + T{2} * fd_delta;
    point_m_1d = y - fd_delta;
    point_m_2d = y - T{2} * fd_delta;

    f_p_1d = function(t, x, point_p_1d, z);
    f_p_2d = function(t, x, point_p_2d, z);
    f_m_1d = function(t, x, point_m_1d, z);
    f_m_2d = function(t, x, point_m_2d, z);
  } else if constexpr (dir == fd_direction::z) {
    point_p_1d = z + fd_delta;
    point_p_2d = z + T{2} * fd_delta;
    point_m_1d = z - fd_delta;
    point_m_2d = z - T{2} * fd_delta;

    f_p_1d = function(t, x, y, point_p_1d);
    f_p_2d = function(t, x, y, point_p_2d);
    f_m_1d = function(t, x, y, point_m_1d);
    f_m_2d = function(t, x, y, point_m_2d);
  }

  return (f_m_2d - T{8} * f_m_1d + T{8} * f_p_1d - f_p_2d) / (T{12} * fd_delta);
}

} // namespace grlensing_tests

#endif // GRLENSING_TESTS_TEST_CONSTANTS_HPP