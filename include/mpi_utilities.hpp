#ifndef GRLENSING_MPI_UTILITIES_HPP
#define GRLENSING_MPI_UTILITIES_HPP

#include <concepts>
#include <limits>
#include <mpi.h>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/iota.hpp>
#include <tuple>
#include <utility>

namespace grlensing {

template <std::unsigned_integral index_t> using index_pair = std::pair<index_t, index_t>;
template <std::unsigned_integral index_t> using index_triplet
    = std::tuple<index_t, index_t, index_t>;

/**
 * Determine the equivalent MPI type of a native data type.
 *
 * @return The MPI::Datatype object representing the underlying data
 */
template <typename T> [[nodiscard]] auto get_mpi_datatype() noexcept(false) -> MPI::Datatype {
  MPI::Datatype mpi_type{MPI::DATATYPE_NULL};

  if constexpr (std::is_same<T, char>::value) {
    mpi_type = MPI::CHAR;
  } else if constexpr (std::is_same<T, signed char>::value) {
    mpi_type = MPI::SIGNED_CHAR;
  } else if constexpr (std::is_same<T, unsigned char>::value) {
    mpi_type = MPI::UNSIGNED_CHAR;
  } else if constexpr (std::is_same<T, wchar_t>::value) {
    mpi_type = MPI::WCHAR;
  } else if constexpr (std::is_same<T, signed short>::value) {
    mpi_type = MPI::SHORT;
  } else if constexpr (std::is_same<T, unsigned short>::value) {
    mpi_type = MPI::UNSIGNED_SHORT;
  } else if constexpr (std::is_same<T, signed int>::value) {
    mpi_type = MPI::INT;
  } else if constexpr (std::is_same<T, unsigned int>::value) {
    mpi_type = MPI::UNSIGNED;
  } else if constexpr (std::is_same<T, signed long int>::value) {
    mpi_type = MPI::LONG;
  } else if constexpr (std::is_same<T, unsigned long int>::value) {
    mpi_type = MPI::UNSIGNED_LONG;
  } else if constexpr (std::is_same<T, signed long long int>::value) {
    mpi_type = MPI::LONG_LONG;
  } else if constexpr (std::is_same<T, unsigned long long int>::value) {
    mpi_type = MPI::UNSIGNED_LONG_LONG;
  } else if constexpr (std::is_same<T, float>::value) {
    mpi_type = MPI::FLOAT;
  } else if constexpr (std::is_same<T, double>::value) {
    mpi_type = MPI::DOUBLE;
  } else if constexpr (std::is_same<T, long double>::value) {
    mpi_type = MPI::LONG_DOUBLE;
  } else if constexpr (std::is_same<T, bool>::value) {
    mpi_type = MPI::BOOL;
  }

  if (mpi_type == MPI::DATATYPE_NULL) {
    std::cout << "From process " << MPI::COMM_WORLD.Get_rank()
              << ": Unable to determine the MPI equivalent of the data type " << typeid(T).name()
              << std::endl;
    throw std::runtime_error("MPI type inference error");
  }

  return mpi_type;
}

/**
 * This function sends data to all eranks except one sender specified rank. This is similar to a
 * broadcast but allows for finer control
 *
 * @param data The data to send.
 * @param sender_rank The rank of the sender process.
 */
template <typename T> auto fine_broadcast(T *data, int sender_rank) noexcept(false) -> void {
  using ranges::views::filter;
  using ranges::views::iota;

  auto exclude_lambda = [&](int i) -> bool { return i != sender_rank; };
  auto proc_list = iota(0, MPI::COMM_WORLD.Get_size()) | filter(exclude_lambda);

  for (int i : proc_list) {
    MPI::COMM_WORLD.Send(static_cast<void *>(data), 1, get_mpi_datatype<T>(), i, 0);
  }
}

template <std::unsigned_integral dest_t, std::signed_integral src_t>
auto safe_cast(src_t src) noexcept(false) -> dest_t {
  if (src < 0 || (std::numeric_limits<dest_t>::max() < std::numeric_limits<src_t>::max())) {
    std::cout << "From process " << MPI::COMM_WORLD.Get_rank() << ": Unable to convert " << src
              << " into the type " << typeid(dest_t).name() << std::endl;
    throw std::runtime_error("Narrowing convertion error");
  } else {
    return dest_t(src);
  }
}

} // namespace grlensing

#endif // GRLENSING_MPI_UTILITIES_HPP