#include "plugin.hpp"

GRLENSING_API grlensing::plugin::plugin(const std::filesystem::path &p) noexcept(false)
    : library(p), get_engine_version_address(nullptr), register_plugin_address(nullptr) {

  /*
   * At this point, the library is guaranteed to be loaded.
   * We can safelly recover the exported functions
   */
  get_engine_version_address
      = library.get_function_pointer<get_engine_version_function>("get_engine_version");
  register_plugin_address
      = library.get_function_pointer<register_plugin_function>("register_plugin");
}