#ifndef GRLENSING_PLUGIN_HPP
#define GRLENSING_PLUGIN_HPP

#include "api_macros.hpp"
#include "shared_library.hpp"

#include <filesystem>

namespace grlensing {

class kernel;

/**
 * Represents a plugin
 */
class plugin {
public:
  /**
   * Load and initialize a plugin
   *
   * @param p Path to the plugin file
   * @return A plugin object.
   */
  GRLENSING_API plugin(const std::filesystem::path &p) noexcept(false);

  /**
   * Queries the plugin for its expected engine version
   *
   * @return A tuple containing the expected version in <major, minor, patch> format
   */
  GRLENSING_API auto get_engine_version() const noexcept -> unsigned {
    return get_engine_version_address();
  }

  /**
   * Register the plugin to a kernel
   *
   * @param kernel Kernel the plugin should register to</param>
   */
  GRLENSING_API auto register_plugin(kernel &kernel) -> void { register_plugin_address(kernel); }

private:
  /**
   * Signature for the version query function
   */
  using get_engine_version_function = unsigned();

  /**
   * Signature for the plugin's registration function
   */
  using register_plugin_function = void(kernel &);

  /**
   * The loaded shared library
   */
  shared_library library;

  /**
   * Function to query for the expected engine version
   */
  get_engine_version_function *get_engine_version_address;

  /**
   * Registers the plugin with the kernel
   */
  register_plugin_function *register_plugin_address;
};

} // namespace grlensing

#endif // GRLENSING_PLUGIN_HPP