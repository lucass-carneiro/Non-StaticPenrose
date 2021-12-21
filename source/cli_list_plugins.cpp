#include "cli.hpp"

void grlensing::list_plugins(grlensing::kernel &kernel) {
  using grlensing::log;
  using grlensing::LogEvent;

  log<LogEvent::message>("Loaded plugins: {:d}", kernel.get_plugin_count());
  kernel.get_storage_server().print_names();
  kernel.get_metric_server().print_names();
}