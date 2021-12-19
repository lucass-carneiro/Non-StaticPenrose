#ifndef GRLENSING_ASCII_WRITER_HPP
#define GRLENSING_ASCII_WRITER_HPP

#include "api_macros.hpp"

#include <cstdio>
#include <string_view>

namespace grlensing {

constexpr const auto std_fclose_wrapper = [](FILE *fp) -> int {
  if (fp != nullptr) {
    std::fclose(fp); // NOLINT
  }
  return EOF;
};

class std_ascii_writer final : public storage_server::trajectory_writer {
public:
  /**
   * Constructs a writer object
   */
  std_ascii_writer()
      : file(nullptr, std_fclose_wrapper), metadata_file(nullptr, std_fclose_wrapper) {}

  /**
   * Tests wether or not the writer can write files in a certain extension type.
   */
  GRLENSING_ASCII_WRITER_API auto can_write_extension(std::string_view extension) -> bool final;

  /**
   * Opens a file for outputing trajectory data.
   *
   * @param p Path to output file
   */
  GRLENSING_ASCII_WRITER_API void open_file(const std::filesystem::path &p,
                                            const std::filesystem::path &m) final;

  /**
   * Opens a file for outputing trajectory data. Does not create a metadata file
   *
   * @param p Path to output file
   * @param m Path to metadata file
   */
  GRLENSING_ASCII_WRITER_API void open_file(const std::filesystem::path &p) final;

  /**
   * Closes an oppened file
   */
  GRLENSING_ASCII_WRITER_API void close_file() final;

  /**
   * Add real value to the next empty column in the data file.
   *
   * @param value The value to be added to the data.
   */
  GRLENSING_ASCII_WRITER_API void push_real(double value) final;

  /**
   * Push the final real value of the row
   *
   * @param value The value to be added to the data.
   */
  GRLENSING_ASCII_WRITER_API void push_final_real(double value) final;

  /**
   * Adds a real valued metadata to the file.
   *
   * @param key The key name of the metadata.
   * @param value The real value of the metadata.
   */
  GRLENSING_ASCII_WRITER_API void push_metadata(std::string_view key, double value) final;

  /**
   *  Adds a unsigned integer valued metadata.
   *
   * @param key The key name of the metadata.
   * @param value The unsigned integer value of the metadata.
   */
  GRLENSING_ASCII_WRITER_API void push_metadata(std::string_view key, std::size_t value) final;

  /**
   * Prints the name of the data writer
   *
   * @return The name of the writer.
   */
  GRLENSING_ASCII_WRITER_API auto name() -> std::string_view final;

private:
  std::unique_ptr<FILE, int (*)(FILE *)> file;
  std::unique_ptr<FILE, int (*)(FILE *)> metadata_file;
};

/**
 * Retrieves the engine version we're expecting
 *
 * @return The engine version the plugin was built against
 */
extern "C" GRLENSING_ASCII_WRITER_API auto get_engine_version() -> unsigned;

/**
 * Registers the plugin to an engine kernel
 * @param kernel Kernel the plugin will be registered to
 */
extern "C" GRLENSING_ASCII_WRITER_API void register_plugin(kernel &kernel);

} // namespace grlensing

#endif // GRLENSING_ASCII_WRITER_HPP