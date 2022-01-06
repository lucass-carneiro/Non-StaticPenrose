#ifndef GRLENSING_MPIMATRIX_HPP
#define GRLENSING_MPIMATRIX_HPP

// TODO: Refactor into map and data

#include "log.hpp"

#include <concepts>
#include <cstdint>
#include <mpi.h>
#include <ranges>
#include <utility>
#include <vector>

namespace grlensing {

enum class ordering_type { row_major, column_major };

template <std::floating_point data_type, std::unsigned_integral index_t = std::uint32_t,
          ordering_type ordering = ordering_type::row_major, bool bound_checked = true>
class mpimatrix {
public:
  using index_pair = std::pair<index_t, index_t>;

  /**
   * Creates a distributed matrix with the given global size.
   *
   * @param rows The global number of rows in the matrix
   * @param cols The global number of columns in the matrix
   */
  mpimatrix(index_t rows, index_t cols)
      : global_rows(rows),
        global_cols(cols),
        global_linear_size(rows * cols),
        owned_liner_range(compute_owned_indexes()),
        owned_linear_size(compute_owned_size()),
        mpi_content_datatype(get_mpi_datatype<data_type>()),
        local_data(owned_linear_size) {}

  /**
   * Determines if the current process owns an index.
   *
   * @param index The index to search.
   * @return The id of the owner process.
   */
  auto is_owner(index_t index) const noexcept(false) -> bool {

    if constexpr (bound_checked) {
      if (index >= global_linear_size) {
        log<LogEvent::error>("The index {} is outside of the matrix.", index);
        throw std::runtime_error("Out of bounds error");
      }
    }

    return (owned_liner_range.first <= index && index <= owned_liner_range.second);
  }

  /**
   * Find which process is the owner of a certain index in the global linear space.
   *
   * This function works by determining if the current process is the owner. If it is, it sends it's
   * rank to all other participating processes. If it is not the owner the function waits from a
   * massage from any other process with the rank of the owner.
   *
   * @param index The index to search for.
   * @return The rank of the owner process.
   */
  auto find_owner(index_t index) const noexcept(false) -> int {
    int owner_rank = -1;

    if (is_owner(index)) {
      owner_rank = MPI::COMM_WORLD.Get_rank();
      send_to_all_but(&owner_rank, owner_rank);
    } else {
      MPI::COMM_WORLD.Recv(&owner_rank, 1, MPI::INT, MPI::ANY_SOURCE, 0);
    }

    if constexpr (bound_checked) {
      if (owner_rank < 0) {
        log<LogEvent::error>("Unable do determine the rank of the owner of index {}.", index);
        throw std::runtime_error("Owner not found");
      }
    }

    return owner_rank;
  }

  /**
   * Converts a global matrix index to global linear index.
   *
   * @param mat_index The global matrix index.
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_pair mat_index) const noexcept(false) -> index_t {
    const auto row = mat_index.first;
    const auto col = mat_index.second;
    index_t I = 0;

    if constexpr (ordering == ordering_type::row_major) {
      I = col + row * global_cols;
    } else if constexpr (ordering == ordering_type::column_major) {
      I = row + col * global_rows;
    }

    if constexpr (bound_checked) {
      if (I >= global_linear_size) {
        log<LogEvent::error>("The index ({},{}) is outside of the matrix.", row, col);
        throw std::runtime_error("Out of bounds error");
      }
    }

    return I;
  }

  /**
   * Converts a global matrix index to global linear index.
   *
   * @param row The global matrix row index.
   * @param column The global matrix row index.
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_t row, index_t col) const noexcept(false) -> index_t {
    return global_matrix_to_linear(std::make_pair(row, col));
  }

  /**
   * Converts a global linear index to global matrix index.
   *
   * @param I The global linear index.
   * @return The corresponding global matrix index pair.
   */
  auto global_linear_to_matrix(index_t I) -> index_pair {

    if constexpr (bound_checked) {
      if (I >= global_linear_size) {
        log<LogEvent::error>("The linear index {} is outside of the matrix.", I);
        throw std::runtime_error("Out of bounds error");
      }
    }

    index_t row = 0, col = 0;

    if constexpr (ordering == ordering_type::row_major) {
      col = I % global_cols;
      row = I / global_cols;
    } else if constexpr (ordering == ordering_type::column_major) {
      row = I % global_rows;
      col = I / global_rows;
    }

    return std::make_pair(row, col);
  }

  /**
   * Gets an element in the matrix.
   *
   * The element is indexed using global indices. If the element does not reside in the current
   * process it is retrieved from the owner process.
   *
   * @param mat_index An index pair in global adressing for the element to be retrieved.
   * @return The value stored at the specified index.
   */
  auto get(index_pair mat_index) const noexcept(false) -> const data_type {
    const auto global_I = global_matrix_to_linear(mat_index);
    data_type retrieved_data = data_type(-1);

    if (is_owner(global_I)) {
      const auto stride = global_I - owned_liner_range.first;
      retrieved_data = local_data.at(stride);

      send_to_all_but(&retrieved_data, MPI::COMM_WORLD.Get_rank());

    } else {
      MPI::COMM_WORLD.Recv(&retrieved_data, 1, mpi_content_datatype, MPI::ANY_SOURCE, 0);
    }

    return retrieved_data;
  }

  /**
   * Gets an element in the matrix.
   *
   * The element is indexed using global indices. If the element does not reside in the current
   * process it is retrieved from the owner process.
   *
   * @param row The row index in global adressing for the element to be retrieved.
   * @param col The column index in global adressing for the element to be retrieved.
   * @return The value stored at the specified index.
   */
  auto get(index_t row, index_t col) const noexcept(false) -> const data_type {
    return get(std::make_pair(row, col));
  }

  /**
   * Sets the data on a given index.
   *
   * @param mat_index The global matrix address of the location to insert.
   * @param value The value to insert.
   */
  void set(index_pair mat_index, data_type value) noexcept(false) {
    const auto global_I = global_matrix_to_linear(mat_index);

    if (is_owner(global_I)) {
      const auto stride = global_I - owned_liner_range.first;
      local_data.at(stride) = value;
    }
  }

  /**
   * Sets the data on a given index.
   *
   * @param row The global matrix row address of the location to insert.
   * @param col The global matrix column address of the location to insert.
   * @param value The value to insert.
   */
  void set(index_t row, index_t col, data_type value) noexcept(false) {
    set(std::make_pair(row, col), value);
  }

  /**
   * TODO: Doc
   */
  auto get_owned_linear_range() const noexcept -> index_pair { return owned_liner_range; }

private:
  using data_vector = std::vector<data_type>;

  template <typename T> auto send_to_all_but(T *data, int excluded_rank) const noexcept(false)
      -> void {
    // Send the retrieved data to all processes but myself
    using std::views::filter;
    using std::views::iota;

    auto exclude_lambda = [&](int i) -> bool { return i != excluded_rank; };
    auto proc_list = iota(0, MPI::COMM_WORLD.Get_size()) | filter(exclude_lambda);

    for (int i : proc_list) {
      MPI::COMM_WORLD.Send(static_cast<void *>(data), 1, get_mpi_datatype<T>(), i, 0);
    }
  }

  /**
   * Compute the pair (start, end) of the locally owned linear indexes
   *
   * These indexes are relative to the global index space in linear ordering.
   *
   * @return the pair (start, end) of the owned indexes.
   */
  auto compute_owned_indexes() const noexcept(false) -> index_pair {
    const auto num_procs = index_t(MPI::COMM_WORLD.Get_size());

    // If there are more processes than arrays, the distribution is ill defined.
    if (num_procs > global_linear_size) {
      log<LogEvent::error>(
          "There are {} MPI processes and {} elements in the distributed matrix. Please make sure "
          "that there are more elements then processes in a distributed matrix.");
      throw std::runtime_error("Load balancing error");
    }

    const auto pid = index_t(MPI::COMM_WORLD.Get_rank());

    /* To make sure that we have an even distribution of points per process, we  first allocate
     * global_linear_size/num_procs for each process. If this division is not exact, the remainder
     * is added to the processes with increasing rank priority. Example:
     *
     * global_index_size = 12
     * num_procs = 5
     * process : 0 1 2 3 4
     * elements: 3 3 2 2 2
     */
    auto base_load = global_linear_size / num_procs;
    const auto load_remainder = global_linear_size % num_procs;

    // The start/end pair is computed.
    auto start_index = index_t(0);

    if (pid < load_remainder) {
      base_load++;
      start_index = pid * base_load;
    } else {
      start_index = pid * base_load + load_remainder;
    }

    const auto end_index = start_index + base_load - 1;
    return std::make_pair(start_index, end_index);
  }

  /**
   * Computes the size of the locally owned portion of the array
   */
  auto compute_owned_size() const noexcept -> index_t {
    return owned_liner_range.second - owned_liner_range.first + 1;
  }

  /**
   * Determine the equivalent MPI type of the underlying data type.
   *
   * @return The MPI::Datatype object representing the underlying data
   */
  template <typename T> [[nodiscard]] auto get_mpi_datatype() const noexcept -> MPI::Datatype {
    if constexpr (std::is_same<T, float>::value) {
      return MPI::FLOAT;
    } else if constexpr (std::is_same<T, double>::value) {
      return MPI::DOUBLE;
    } else {
      return MPI::LONG_DOUBLE;
    }
  }

  /**
   * The global numbr of rows in the matrix.
   */
  index_t global_rows;

  /**
   * The global numbr of columns in the matrix.
   */
  index_t global_cols;

  /**
   * The total size of the matrix in linear coordinates
   */
  index_t global_linear_size;

  /**
   * The (start,end) index par of the locally owned portion of the array, given in linear
   * coordinates
   */
  index_pair owned_liner_range;

  /**
   * The number of linear indices owned by this portion of the array.
   */
  index_t owned_linear_size;

  /**
   * The MPI data type associated with the underlying data.
   */
  const MPI::Datatype mpi_content_datatype;

  /**
   * Linear array with the local data
   */
  data_vector local_data;
};

} // namespace grlensing

#endif // GRLENSING_MPIMATRIX_HPP