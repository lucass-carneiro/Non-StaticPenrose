#ifndef GRLENSING_TENSOR_TYPES_HPP
#define GRLENSING_TENSOR_TYPES_HPP

#include <array>
#include <cstddef>

namespace grlensing {

template <std::size_t d1 = 3, typename T = double> struct callable_array {
  std::array<T, d1> data;

  [[nodiscard]] inline constexpr auto operator()(std::size_t i) const -> T { return data[i]; }

  inline auto operator[](std::size_t i) -> T & { return data[i]; }
  inline auto operator[](std::size_t i) const -> const T & { return data[i]; }
};

template <std::size_t d1, typename T> auto operator*(T &&a, const callable_array<d1, T> &b)
    -> callable_array<d1, T> {
  callable_array<d1, T> ret{};

  for (std::size_t i = 0; i < d1; i++) {
    ret[i] = a * b(i);
  }

  return ret;
}

template <std::size_t d1, typename T>
auto operator+(const callable_array<d1, T> &a, const callable_array<d1, T> &b)
    -> const callable_array<d1, T> {
  callable_array<d1, T> ret{};

  for (std::size_t i = 0; i < d1; i++) {
    ret[i] = a(i) + b(i);
  }

  return ret;
}

template <std::size_t d1 = 3, std::size_t d2 = 3, typename T = double> struct callable_matrix {
  callable_array<d1, callable_array<d2, T>> data;

  [[nodiscard]] inline constexpr auto operator()(std::size_t i, std::size_t j) const -> T {
    return data[i][j];
  }

  inline auto operator[](std::size_t i) -> callable_array<d2, T> & { return data[i]; }
  inline auto operator[](std::size_t i) const -> const callable_array<d2, T> & { return data[i]; }
};

template <std::size_t d1, std::size_t d2, typename T>
auto operator*(T &&a, const callable_matrix<d1, d2, T> &b) -> callable_matrix<d1, d2, T> {
  callable_matrix<d1, d2, T> ret{};

  for (std::size_t i = 0; i < d1; i++) {
    for (std::size_t j = 0; j < d2; j++) {
      ret[i][j] = a * b(i, j);
    }
  }

  return ret;
}

template <std::size_t d1, std::size_t d2, typename T>
auto operator+(const callable_matrix<d1, d2, T> &a, const callable_matrix<d1, d2, T> &b)
    -> const callable_matrix<d1, d2, T> {
  callable_matrix<d1, d2, T> ret{};

  for (std::size_t i = 0; i < d1; i++) {
    for (std::size_t j = 0; j < d2; j++) {
      ret[i][j] = a(i, j) + b(i, j);
    }
  }

  return ret;
}

template <std::size_t d1 = 3, std::size_t d2 = 3, std::size_t d3 = 3, typename T = double>
struct callable_3tensor {
  std::array<callable_matrix<d2, d3, T>, d1> data;

  [[nodiscard]] inline constexpr auto operator()(std::size_t i, std::size_t j, std::size_t k) const
      -> T {
    return data[i](j, k);
  }

  inline auto operator[](std::size_t i) -> callable_matrix<d2, d3, T> & { return data[i]; }
  inline auto operator[](std::size_t i) const -> const callable_matrix<d2, d3, T> & {
    return data[i];
  }
};

} // namespace grlensing

#endif // GRLENSING_TENSOR_TYPES_HPP