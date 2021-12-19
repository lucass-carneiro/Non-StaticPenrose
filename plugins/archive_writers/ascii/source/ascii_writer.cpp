#include "ascii_writer.hpp"

// clang-format off
#include <cstdio>
#include <fmt/core.h>

#include <cerrno>
#include <cstring>
#include <stdexcept>

using grlensing::kernel;
using grlensing::storage_server;
using grlensing::std_ascii_writer;
// clang-format on

GRLENSING_ASCII_WRITER_API auto std_ascii_writer::can_write_extension(std::string_view extension)
    -> bool {
  if (extension == ".ascii" || extension == ".dat" || extension == ".txt")
    return true;
  else
    return false;
}

GRLENSING_ASCII_WRITER_API void
std_ascii_writer::open_file(const std::filesystem::path &p,
                            const std::filesystem::path &m) noexcept(false) {

  /* Precondition: Requires that the pointers are null prior to oppening. If this is false, it means
   * that we are trying to open a new file without closing the ones that are already open. This
   * would cause a memory leak. Open multiple files at once is likelly a bug im the user code, an
   * exception is thrown to halt the code. Here, file and metadata_file are not expected to change
   * independently, so we could chek either one of them. For redundancy, both pointers are tested.
   */
  if (file != nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing data because a file is already "
                         "open. Close all files before opening a new one",
                         p.string());
    throw std::runtime_error("IO error");
  }

  if (metadata_file != nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing data because a file is already "
                         "open. Close all files before opening a new one",
                         m.string());
    throw std::runtime_error("IO error");
  }

  // Open files
  file = std::unique_ptr<FILE, int (*)(FILE *)>(std::fopen(p.c_str(), "w"), std_fclose_wrapper);
  metadata_file
      = std::unique_ptr<FILE, int (*)(FILE *)>(std::fopen(m.c_str(), "w"), std_fclose_wrapper);

  // Postcondition: Ensures that the fopen operation succedded.
  if (file == nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing data: ", p.string(),
                         std::strerror(errno));
    throw std::runtime_error("IO error");
  }

  if (metadata_file == nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing metadata: ", m.string(),
                         std::strerror(errno));
    throw std::runtime_error("IO error");
  }
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::open_file(const std::filesystem::path &p) {
  /* Precondition: Requires that the pointers are null prior to oppening. If this is false, it means
   * that we are trying to open a new file without closing the ones that are already open. This
   * would cause a memory leak. Open multiple files at once is likelly a bug im the user code, an
   * exception is thrown to halt the code. Here, file and metadata_file are not expected to change
   * independently, so we could chek either one of them. For redundancy, both pointers are tested.
   */
  if (file != nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing data because a file is already "
                         "open. Close all files before opening a new one",
                         p.string());
    throw std::runtime_error("IO error");
  }

  // Open files
  file = std::unique_ptr<FILE, int (*)(FILE *)>(std::fopen(p.c_str(), "w"), std_fclose_wrapper);

  // Postcondition: Ensures that the fopen operation succedded.
  if (file == nullptr) {
    log<LogEvent::error>("Unable to open file {:s} for outputing data: ", p.string(),
                         std::strerror(errno));
    throw std::runtime_error("IO error");
  }
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::close_file() {
  file.reset(nullptr);
  metadata_file.reset(nullptr);
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::push_real(double value) {
  /* Precondition: Ensure that the file is open before writing to it. Writing to a closed file is
   * likelly a bug in the user code so an exception is thrown
   */
  if (file == nullptr) {
    log<LogEvent::error>("Unable to write data to a closed file stream.");
    throw std::runtime_error("IO error");
  } else {
    fmt::print(file.get(), "{0:.16e}    ", value);
  }
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::push_final_real(double value) {
  /* Precondition: Ensure that the file is open before writing to it. Writing to a closed file is
   * likelly a bug in the user code so an exception is thrown
   */
  if (file == nullptr) {
    log<LogEvent::error>("Unable to write data to a closed file stream.");
    throw std::runtime_error("IO error");
  } else {
    fmt::print(file.get(), "{0:.16e}\n", value);
  }
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::push_metadata(std::string_view key,
                                                                double value) {
  /* Precondition: Ensure that the file is open before writing to it. Writing to a closed file is
   * likelly a bug in the user code so an exception is thrown
   */
  if (metadata_file == nullptr) {
    log<LogEvent::error>("Unable to write data to a closed file stream.");
    throw std::runtime_error("IO error");
  } else {
    fmt::print(metadata_file.get(), "{:s} = ", key);
    fmt::print(metadata_file.get(), "{0:.16e}\n", value);
  }
}

GRLENSING_ASCII_WRITER_API void std_ascii_writer::push_metadata(std::string_view key,
                                                                std::size_t value) {
  /* Precondition: Ensure that the file is open before writing to it. Writing to a closed file is
   * likelly a bug in the user code so an exception is thrown
   */
  if (metadata_file == nullptr) {
    log<LogEvent::error>("Unable to write data to a closed file stream.");
    throw std::runtime_error("IO error");
  } else {
    fmt::print(metadata_file.get(), "{:s} = ", key);
    fmt::print(metadata_file.get(), "{:d}\n", value);
  }
}

GRLENSING_ASCII_WRITER_API auto std_ascii_writer::name() -> std::string_view {
  return "Default trajectory ASCII file writer.";
}

extern "C" GRLENSING_ASCII_WRITER_API auto get_engine_version() -> unsigned { return unsigned(1); }

extern "C" GRLENSING_ASCII_WRITER_API void register_plugin(kernel &kernel) {
  kernel.get_storage_server().add_trajectory_writer(
      std::unique_ptr<storage_server::trajectory_writer>(new std_ascii_writer()));
}