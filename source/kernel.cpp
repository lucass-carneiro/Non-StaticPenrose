#include "kernel.hpp"

#include "log.hpp"

GRLENSING_API void grlensing::kernel::load_plugin(const std::filesystem::path &p) {
  if (loaded_plugins.find(p.filename().string()) == loaded_plugins.end()) {
    const auto result = loaded_plugins.insert({p.filename().string(), plugin(p)});
    if (result.second) {
      result.first->second.register_plugin(*this);
    } else {
      log<LogEvent::error>("Unable to insert the plugin {:s} in the plugin registry.",
                           p.filename().string());
      throw std::runtime_error("Plugin loading error");
    }
  }
}