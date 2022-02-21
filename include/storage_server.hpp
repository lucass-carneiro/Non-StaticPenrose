#ifndef GRLENSING_STORAGE_SERVER_HPP
#define GRLENSING_STORAGE_SERVER_HPP

#include "api_macros.hpp"
#include "log.hpp"

#include <filesystem>
#include <list>
#include <string_view>

namespace grlensing {

// ----------------------------------------------------------------------- //

class archive {}; // Dummy

// ----------------------------------------------------------------------- //

/**
 * Manages storage related stuff
 */
class storage_server {
public:
  /**
   * Reads archive files like zips and rars.
   * Abstract interface to be implemented by the user
   */
  class archive_reader { // TODO: add remaining constructors
  public:
    /**
     * Checks whether the reader can open an archive
     *
     * @param p Path of the file that will be checked
     */
    virtual auto can_open_archive(const std::filesystem::path &p) -> bool = 0;

    /**
     * Opens an archive for reading</summary>
     *
     * @param p Path of the file to open
     */
    virtual auto open_archive(const std::filesystem::path &p) -> std::unique_ptr<archive> = 0;

    /**
     * The display name of the archive reader
     *
     * @return A string with the name of the reader.
     */
    virtual auto name() -> std::string_view = 0;

    virtual ~archive_reader() = default;
  };

  /**
   * Writes files containing the result of trajectory integration.
   * Abstract interface to be implemented by the user.
   */
  class trajectory_writer { // TODO: add remaining constructors
  public:
    /**
     * Tests wether or not the writer can write files in a certain extension type.
     */
    virtual auto can_write_extension(std::string_view extension) -> bool = 0;

    /**
     * Opens a file for outputing trajectory data.
     *
     * @param p Path to output file
     * @param m Path to metadata file
     */
    virtual void open_file(const std::filesystem::path &p, const std::filesystem::path &m) = 0;

    /**
     * Opens a file for outputing trajectory data. Does not create a metadata file
     *
     * @param p Path to output file
     */
    virtual void open_file(const std::filesystem::path &p) = 0;

    /**
     * Closes oppened files.
     */
    virtual void close_file() = 0;

    /**
     * Add real value to the next empty column in the data file.
     *
     * @param value The value to be added to the data.
     */
    virtual void push_real(double value) = 0;

    /**
     * Push the final real value of the row
     *
     * @param value The value to be added to the data.
     */
    virtual void push_final_real(double value) = 0;

    /**
     * Adds a real valued metadata to the file.
     *
     * @param key The key name of the metadata.
     * @param value The real value of the metadata.
     */
    virtual void push_metadata(std::string_view key, double value) = 0;

    /**
     *  Adds a unsigned integer valued metadata.
     *
     * @param key The key name of the metadata.
     * @param value The unsigned integer value of the metadata.
     */
    virtual void push_metadata(std::string_view key, std::size_t value) = 0;

    /**
     * Prints the name of the data writer
     *
     * @return The name of the writer.
     */
    virtual auto name() -> std::string_view = 0;

    virtual ~trajectory_writer() = default;
  };

  /**
   * Allows plugins to add new archive readers
   *
   * @param reader Pointer to a archive reader.
   */
  GRLENSING_API void add_archive_reader(std::unique_ptr<archive_reader> reader) {
    archive_readers.push_back(std::move(reader));
  }

  /**
   * Allows plugins to add new trajectory writers.
   *
   * @param writer Pointer to file writer
   */
  GRLENSING_API void add_trajectory_writer(std::unique_ptr<trajectory_writer> writer) {
    writers.push_back(std::move(writer));
  }

  /**
   * Opens an archive by searching for a matching archive reader
   *
   * @param p Path to a file to be oppened
   */
  GRLENSING_API auto open_archive(const std::filesystem::path &p) -> std::unique_ptr<archive>;

  /**
   * Gets a writer capable of handling a certain file extension
   *
   * @param extension The extension that the writer should handle.
   * @return A constant reference to a unique_ptr wraping the writer.
   */
  GRLENSING_API auto get_writer(std::string_view extension) const noexcept(false)
      -> const std::unique_ptr<trajectory_writer> &;

  /**
   * Print the name of the storage servers, as provided by the plugin
   * implementer.
   */
  GRLENSING_API void print_names() const;

private:
  /**
   * A list of archive readers
   */
  using archive_reader_list = std::list<std::unique_ptr<archive_reader>>;

  /**
   * All available archive readers
   */
  archive_reader_list archive_readers;

  /**
   * A list of tajectory writers
   */
  using writer_list = std::list<std::unique_ptr<trajectory_writer>>;

  /**
   * All available tajectory writers
   */
  writer_list writers;
};

} // namespace grlensing

#endif // GRLENSING_STORAGE_SERVER_HPP