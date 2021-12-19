#ifndef GRLENSING_KERNEL_HPP
#define GRLENSING_KERNEL_HPP

#include "api_macros.hpp"
#include "metric_server.hpp"
#include "plugin.hpp"
#include "storage_server.hpp"

#include <filesystem>
#include <map>
#include <string>

namespace grlensing {

/**
 * The engine core and plugin manager
 */
class kernel {
public:
  /**
   * Accesses the storage server
   */
  GRLENSING_API auto get_storage_server() -> storage_server & { return stg_serv; }

  /**
   * Accesses the metric server
   */
  GRLENSING_API auto get_metric_server() -> metric_server & { return mtr_serv; }

  /**
   * Loads a plugin
   *
   * @param p Path to a plugin file
   */
  GRLENSING_API void load_plugin(const std::filesystem::path &p);

  /**
   * Returns the plugin count
   *
   * @return The plugin count.
   */
  GRLENSING_API auto get_plugin_count() -> std::size_t { return loaded_plugins.size(); }

private:
  /**
   * Map of plugins by their associated file names
   */
  using plugin_map = std::map<std::string, plugin>;

  /**
   * All plugins currently loaded
   */
  plugin_map loaded_plugins;

  /**
   * Manages storage-related tasks for the engine
   */
  storage_server stg_serv;

  /**
   * Manages all the loaded metrics
   */
  metric_server mtr_serv;
};

} // namespace grlensing

#endif // GRLENSING_KERNEL_HPP