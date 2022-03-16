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

template <std::size_t d1 = 3, std::size_t d2 = 3, typename T = double> struct callable_matrix {
  std::array<std::array<T, d2>, d1> data;

  [[nodiscard]] inline constexpr auto operator()(std::size_t i, std::size_t j) const -> T {
    return data[i][j];
  }

  inline auto operator[](std::size_t i) -> std::array<T, d2> & { return data[i]; }
  inline auto operator[](std::size_t i) const -> const std::array<T, d2> & { return data[i]; }
};

template <std::size_t d1 = 3, std::size_t d2 = 3, std::size_t d3 = 3, typename T = double>
struct callable_3tensor {
  std::array<std::array<std::array<T, d3>, d2>, d1> data;

  [[nodiscard]] inline constexpr auto operator()(std::size_t i, std::size_t j, std::size_t k) const
      -> T {
    return data[i][j][k];
  }

  inline auto operator[](std::size_t i) -> std::array<std::array<T, d3>, d2> & { return data[i]; }
  inline auto operator[](std::size_t i) const -> const std::array<std::array<T, d3>, d2> & {
    return data[i];
  }
};

} // namespace grlensing

#endif // GRLENSING_TENSOR_TYPES_HPP