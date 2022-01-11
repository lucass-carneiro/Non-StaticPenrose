#ifndef GRLENSING_MPI_INDEX_MAP_3D_HPP
#define GRLENSING_MPI_INDEX_MAP_3D_HPP

#include "mpi_utilities.hpp"

#include <cstdint>

namespace grlensing {

enum class ordering_type { row_major, column_major };

/**
 * Creates a index map for distributed computations
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
 * $(i_1,...,i_{d-1})$ is reduced to a  single value $I$ that is uniquelly invertible back to the
 * origin tuple. The adress space of $I$ is called the linear address space.
 *
 * 4. Data adress space: The space of the indexs in the tuple $(i_1,...,i_{d-1})$
 */
template <std::unsigned_integral index_t = std::uint32_t,
          ordering_type ordering = ordering_type::row_major, bool bound_checked = true>
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
      : global_data_Nx(Nx),
        global_data_Ny(Ny),
        global_data_Nz(Nz),
        global_size(global_data_Nx * global_data_Ny * global_data_Nz),
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
   * Converts a global data space index to a global linear index.
   *
   * @param data_index The global matrix index.
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_triplet<index_t> data_index) const noexcept(false) -> index_t {
    const auto i = std::get<0>(data_index);
    const auto j = std::get<1>(data_index);
    const auto k = std::get<2>(data_index);
    auto I = index_t{0};

    if constexpr (ordering == ordering_type::row_major) {
      I = k + (j + i * global_data_Ny) * global_data_Nz;
    } else if constexpr (ordering == ordering_type::column_major) {
      I = i + (j + k * global_data_Ny) * global_data_Nx;
    }

    if constexpr (bound_checked) {
      if (I >= global_size) {
        std::cout << "From process " << MPI::COMM_WORLD.Get_rank() << ": The index (" << i << ","
                  << j << "," << k << ") is out of bounds" << std::endl;
        throw std::runtime_error("Out of bounds error");
      }
    }

    return I;
  }

  /**
   * Converts a global data space index to a global linear index.
   *
   * @param data_index The global matrix index.
   * @return The corresponding global linear index.
   */
  auto global_matrix_to_linear(index_t i, index_t j, index_t k) const noexcept(false) -> index_t {
    return global_matrix_to_linear(std::make_tuple(i, j, k));
  }

  /**
   * Converts a global linear index to global matrix index.
   *
   * @param I The global linear index.
   * @return The corresponding global index triplet.
   */
  auto global_linear_to_matrix(index_t I) -> index_triplet<index_t> {
    if constexpr (bound_checked) {
      if (I >= global_size) {
        std::cout << "From process " << MPI::COMM_WORLD.Get_rank() << ": The index " << I
                  << " is out of bounds" << std::endl;
        throw std::runtime_error("Out of bounds error");
      }
    }

    index_t i = 0, j = 0, k = 0;

    if constexpr (ordering == ordering_type::row_major) {
      k = I % global_data_Nz;
      I /= global_data_Nz;

      j = I % global_data_Ny;
      I /= global_data_Ny;

      i = I;
    } else if constexpr (ordering == ordering_type::column_major) {
      i = I % global_data_Nx;
      I /= global_data_Nx;

      j = I % global_data_Ny;
      I /= global_data_Ny;

      k = I;
    }

    return std::make_tuple(i, j, k);
  }

  /**
   * Determines if the current process owns a global linear index.
   *
   * @param gli The global linear index to search.
   * @return A boolean indicating if the current process is the owner of the index.
   */
  auto is_owner(index_t gli) const noexcept(false) -> bool {

    if constexpr (bound_checked) {
      if (gli >= global_size) {
        std::cout << "From process " << MPI::COMM_WORLD.Get_rank() << ": The index (" << gli
                  << " is out of bounds" << std::endl;
        throw std::runtime_error("Out of bounds error");
      }
    }

    return (global_linear_owned_range.first <= gli && gli <= global_linear_owned_range.second);
  }

  /**
   * Find which process is the owner of a certain index given in the global linear coordinates.
   *
   * This function works by first determining if the current process is the owner.
   * If it is, it sends it's rank to all other participating processes. If it is not it waits for a
   * message from any other process. A messege is guaranteed to be recieved because either some
   * process is the owner or an exception is thrown.
   *
   * @param index The index to search for.
   * @return The rank of the owner process.
   */
  auto find_owner(index_t index) const noexcept(false) -> int {
    int owner_rank = -1;

    if (is_owner(index)) {
      owner_rank = MPI::COMM_WORLD.Get_rank();
      fine_broadcast(&owner_rank, owner_rank);
    } else {
      MPI::COMM_WORLD.Recv(&owner_rank, 1, MPI::INT, MPI::ANY_SOURCE, 0);
    }

    if constexpr (bound_checked) {
      if (owner_rank < 0) {
        std::cout << "From process " << MPI::COMM_WORLD.Get_rank()
                  << "Unable do determine the rank of the owner of index " << index << std::endl;
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
   * Computes the global_linear_owned_range. See the variable's documentation for further details
   *
   * @return The pair [start, end] of the owned indexes.
   */
  auto compute_glor() const noexcept(false) -> index_pair<index_t> {
    const auto num_procs = safe_cast<index_t>(MPI::COMM_WORLD.Get_size());
    const auto rank = safe_cast<index_t>(MPI::COMM_WORLD.Get_rank());

    if (num_procs > global_size) {
      std::cout << "From process " << MPI::COMM_WORLD.Get_rank() << ": There are " << num_procs
                << " MPI processes and " << global_size
                << "elements in the data structure. Please make sure that there are more elements "
                   "then processes."
                << std::endl;
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
    auto base_load = global_size / num_procs;
    const auto load_remainder = global_size % num_procs;

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
  index_t global_size;

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