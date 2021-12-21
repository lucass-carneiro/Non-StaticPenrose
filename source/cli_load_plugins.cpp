#include "cli.hpp"

void grlensing::push_to_paths(const std::filesystem::path &plugin_folder, const YAML::Node &array,
                              const std::string &suffix,
                              std::vector<std::filesystem::path> &path_vector) {
  for (const auto &plugin_name : array) {
    const auto base_name = plugin_name.as<std::string>();

    auto lib_name = std::string("libgrlensing_");
    lib_name += base_name;
    lib_name += suffix;

    const auto lib_path = std::filesystem::path(lib_name);
    path_vector.push_back(plugin_folder / lib_path);
  }
}

auto grlensing::make_plugin_paths(const YAML::Node &config_file)
    -> std::vector<std::filesystem::path> {
  // Base plugin folders
  const auto metric_plugin_folder = std::filesystem::path(
      config_file["plugin_settings"]["metric_plugin_folder"].as<std::string>());
  const auto archive_reader_plugin_folder = std::filesystem::path(
      config_file["plugin_settings"]["archive_reader_plugin_folder"].as<std::string>());
  const auto archive_writer_plugin_folder = std::filesystem::path(
      config_file["plugin_settings"]["archive_writer_plugin_folder"].as<std::string>());

  // Find all plugins in the configured paths
  std::vector<std::filesystem::path> plugin_paths;
  plugin_paths.reserve(config_file["plugin_settings"]["load_metrics"].size()
                       + config_file["plugin_settings"]["load_readers"].size());

  push_to_paths(metric_plugin_folder, config_file["plugin_settings"]["load_metrics"], "_metric.so",
                plugin_paths);
  push_to_paths(archive_reader_plugin_folder, config_file["plugin_settings"]["load_readers"],
                "_archive.so", plugin_paths);
  push_to_paths(archive_writer_plugin_folder, config_file["plugin_settings"]["load_writers"],
                "_writer.so", plugin_paths);

  return plugin_paths;
}

void grlensing::load_plugins(grlensing::kernel &kernel, const YAML::Node &config_file) {
  const auto plugin_paths = make_plugin_paths(config_file);

  for (const auto &plugin_path : plugin_paths) {
    kernel.load_plugin(plugin_path);
  }
}