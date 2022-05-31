#ifndef GRLENSING_MPI_INDEX_MAP_3D_HPP
#define GRLENSING_MPI_INDEX_MAP_3D_HPP

#include "log.hpp"
#include "mpi_utilities.hpp"

#include <cstdint>

namespace grlensing {

enum class ordering_type { row_major, column_major };

/**
 * Creates a 3D index map for distributed computations
 *
 * Nomenclature:
 *
 * 1. Global index: An index that lables something globally, independent of the data distribution
 * layout. If a communicator has a single process, the local and global indices coincide.
 *
 * 2. Local index: An index that labels something on the current process. Local and global indices
 * do not coincide unless there is only one process in the communicatore
 *
 * 3. Linear adress space: If a data structure has $d$ dimentions, the index tuple
 * $(i_0,...,i_{d-1})$ is reduced to a  single value $I$ that is uniquelly invertible back to the
 * origin tuple. The adress space of $I$ is called the linear address space.
 *
 * 4. Data adress space: The space of the indexs in the tuple $(i_0,...,i_{d-1})$
 *
 * @tparam index_t The type to be used in index computations
 * @tparam ordering The ordering to use when flattening indexes (row major vs column major)
 * @tparam safe_operations Whether or not to perform oerflow and bound checking.
 */
template <std::unsigned_integral index_t = std::uint32_t,
          ordering_type ordering = ordering_type::row_major, bool safe_operations = true>
class mpi_index_map_3D {
public:
  /**
   * Creates a distributed index map with the given 3 dimentional structure.
   *
   * @param Nx The global number of elements in the x direction
   * @param Ny The global number of elements in the y direction
   * @param Nz The global number of elements in the z direction
   */
  mpi_index_map_3D(index_t Nx, index_t Ny, index_t Nz)
      : global_data_Nx{Nx},
        global_data_Ny{Ny},
        global_data_Nz{Nz},
        global_linear_size(compute_global_size()),
        global_linear_owned_range(compute_glor()),
        global_linear_owned_size(compute_glos()) {}

  /**
   * Creates a distributed index map with the given 3 dimentional structure.
   *
   * @param triplet A triplet with the global number of indices in each direction
   */
  mpi_index_map_3D(index_triplet<index_t> triplet)
      : mpi_index_map_3D(std::get<0>(triplet), std::get<1>(triplet), std::get<2>(triplet)) {}

  /**
   * Converts a global data space index to a global linear space index.
   *
   * @param data_index The global matrix index.
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_triplet<index_t> data_index) const noexcept(false) -> index_t {
    const auto i = std::get<0>(data_index);
    const auto j = std::get<1>(data_index);
    const auto k = std::get<2>(data_index);
    index_t I{0};

    if constexpr (ordering == ordering_type::column_major) {
      I = k + global_data_Nz * (j + global_data_Ny * i);
    } else if constexpr (ordering == ordering_type::row_major) {
      I = i + global_data_Nx * (j + global_data_Ny * k);
    }

    if constexpr (safe_operations) {
      if (I >= global_linear_size) {
        log<LogEvent::error>("The index ({},{},{}) is out of bounds", i, j, k);
        throw std::runtime_error("Out of bounds error");
      }
    }

    return I;
  }

  /**
   * Converts a global data space index to a global linear index.
   *
   * @param data_index The global matrix index to "flatten".
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_t i, index_t j, index_t k) const noexcept(false) -> index_t {
    return global_matrix_to_linear(std::make_tuple(i, j, k));
  }

  /**
   * Converts a global linear index to global matrix index.
   *
   * @param gli The global linear index.
   * @return The corresponding global index triplet.
   */
  auto global_linear_to_global_matrix(index_t gli) const -> index_triplet<index_t> {

    if constexpr (safe_operations) {
      if (gli >= global_linear_size) {
        log<LogEvent::error>("The index {} is out of bounds", gli);
        throw std::runtime_error("Out of bounds error");
      }
    }

    index_t i = 0, j = 0, k = 0;

    if constexpr (ordering == ordering_type::column_major) {
      k = gli % global_data_Nz;
      gli /= global_data_Nz;

      j = gli % global_data_Ny;
      gli /= global_data_Ny;

      i = gli;
    } else if constexpr (ordering == ordering_type::row_major) {
      i = gli % global_data_Nx;
      gli /= global_data_Nx;

      j = gli % global_data_Ny;
      gli /= global_data_Ny;

      k = gli;
    }

    return std::make_tuple(i, j, k);
  }

  /**
   * Determines if the current process owns a global linear index.
   *
   * @param gli The global linear index to search in the current process.
   * @return A boolean indicating if the current process is the owner of the index.
   */
  auto is_owner(index_t gli) const noexcept(false) -> bool {

    if constexpr (safe_operations) {
      if (gli >= global_linear_size) {
        log<LogEvent::error>("The index {} is out of bounds", gli);
        throw std::runtime_error("Out of bounds error");
      }
    }

    return (global_linear_owned_range.first <= gli && gli <= global_linear_owned_range.second);
  }

  /**
   * Find which process is the owner of a certain global linear index.
   *
   * This function works by first determining if the current process is the owner.
   * If it is, it sends it's rank to all other participating processes. If it is not it waits for a
   * message from any other process.
   *
   * @param gli The global linear index to search for.
   * @return The rank of the owner process.
   */
  auto find_owner(index_t gli) const noexcept(false) -> int {
    int owner_rank = -1;

    if (is_owner(gli)) {
      owner_rank = MPI::COMM_WORLD.Get_rank();
      fine_broadcast(&owner_rank, owner_rank);
    } else {
      MPI::COMM_WORLD.Recv(&owner_rank, 1, MPI::INT, MPI::ANY_SOURCE, 0);
    }

    if constexpr (safe_operations) {
      if (owner_rank < 0) {
        log<LogEvent::error>("Unable do determine the rank of the owner of index {}", gli);
        throw std::runtime_error("Owner not found");
      }
    }

    return owner_rank;
  }

  /**
   * Returns the owned linear range in global coordinates
   *
   * @return The owned linear range in global coordinates
   */
  auto get_glor() const noexcept -> index_pair<index_t> { return global_linear_owned_range; }

private:
  /**
   * Computes total the amount of linear indexes
   *
   * If safe_operations is set to true, this function checks for overflow while performing the
   * necessary multiplication using compiler builtins.
   */
  auto compute_global_size() -> index_t {
    if constexpr (safe_operations) {
      index_t Nx_times_Ny{0};
      index_t Nz_times_Nx_times_Ny{0};

      if (__builtin_mul_overflow(global_data_Nx, global_data_Ny, &Nx_times_Ny)) {
        log<LogEvent::error>("Overflow occurrend while computing Nx and Ny");
        throw std::runtime_error("Overflow error");
      }

      if (__builtin_mul_overflow(global_data_Nz, Nx_times_Ny, &Nz_times_Nx_times_Ny)) {
        log<LogEvent::error>("Overflow occurrend while computing Nz * (Nx * Ny)");
        throw std::runtime_error("Overflow error");
      }
    }

    return global_data_Nx * global_data_Ny * global_data_Nz;
  }

  /**
   * Computes the global_linear_owned_range. See the variable's documentation for further details
   *
   * @return The pair [start, end] of the owned indexes.
   */
  auto compute_glor() const noexcept(false) -> index_pair<index_t> {
    const auto num_procs = safe_cast<index_t>(MPI::COMM_WORLD.Get_size());
    const auto rank = safe_cast<index_t>(MPI::COMM_WORLD.Get_rank());

    if (num_procs > global_linear_size) {
      log<LogEvent::error>("There are {} MPI processes and {} elements in the data structure. Make "
                           "sure that there are more elements than processes.",
                           num_procs, global_linear_size);
      throw std::runtime_error("Load balancing error");
    }

    /* To make sure that we have an even distribution of points per process, we first assign
     * global_size/num_procs points to each process. If this division is not exact, the remainder
     * is added to all process linearly starting with rank 0. Example:
     *
     * global_size = 12
     * num_procs = 5
     * rank:     0 1 2 3 4
     * elements: 3 3 2 2 2
     */
    auto base_load = global_linear_size / num_procs;
    const auto load_remainder = global_linear_size % num_procs;

    auto start_index = index_t{0};

    if (rank < load_remainder) {
      base_load++;
      start_index = rank * base_load;
    } else {
      start_index = rank * base_load + load_remainder;
    }

    const auto end_index = start_index + base_load - 1;
    return std::make_pair(start_index, end_index);
  }

  /**
   * Computes global_linear_owned_size. See the variables documentaion for further information.
   *
   * @return The global linear owned size of the data structure.
   */
  auto compute_glos() const noexcept -> index_t {
    return global_linear_owned_range.second - global_linear_owned_range.first + 1;
  }

  /**
   * The global number of elements in the "x direction" of the data structure.
   */
  index_t global_data_Nx;

  /**
   * The global number of elements in the "x direction" of the data structure.
   */
  index_t global_data_Ny;

  /**
   * The global number of elements in the "x direction" of the data structure.
   */
  index_t global_data_Nz;

  /**
   * The total size of the data structure in linear
   */
  index_t global_linear_size;

  /**
   * The range of global indices, in linear adressing, owned by a process, as a (inclusive at the
   * end values) pair in the form [start, end]
   */
  index_pair<index_t> global_linear_owned_range;

  /**
   * The number of linear indices owned by this portion of the array.
   */
  index_t global_linear_owned_size;
};

} // namespace grlensing

#endif // GRLENSING_MPI_INDEX_MAP_3D_HPP