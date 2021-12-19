
//#define GRLENSING_ZIP_PLUGIN_SOURCE

#include "api_macros.hpp"

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace grlensing {

/**
 * Zip archive reader
 */
class zip_archive_reader : public storage_server::archive_reader {
  /**
   * Checks whether the archive reader can open an archive
   *
   * @param p Archive that will be tested</param>
   */
public:
  GRLENSING_ZIP_PLUGIN_API zip_archive_reader() = default;

  GRLENSING_ZIP_PLUGIN_API bool can_open_archive(const std::filesystem::path &p) {
    return p.extension() == ".zip";
  }

  /**
   * Opens an archive for reading
   * @param p Archive that will be opened
   */
  GRLENSING_ZIP_PLUGIN_API std::unique_ptr<archive> open_archive(const std::filesystem::path &p) {
    if (!can_open_archive(p)) {
      throw std::runtime_error("No Zip archive");
    }

    return std::unique_ptr<archive>(new archive());
  }

  GRLENSING_ZIP_PLUGIN_API std::string_view name() { return "Default ZIP archive handler"; }
};

/**
 * Retrieves the engine version we're expecting
 *
 * @return The engine version the plugin was built against
 */
extern "C" GRLENSING_ZIP_PLUGIN_API unsigned get_engine_version() { return unsigned(1); }

/**
 * Registers the plugin to an engine kernel
 * @param kernel Kernel the plugin will be registered to
 */
extern "C" GRLENSING_ZIP_PLUGIN_API void register_plugin(kernel &kernel) {
  kernel.get_storage_server().add_archive_reader(
      std::unique_ptr<storage_server::archive_reader>(new zip_archive_reader()));
}

} // namespace grlensing