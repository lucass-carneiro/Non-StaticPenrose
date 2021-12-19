#include "shared_library.hpp"

auto grlensing::shared_library::validate_libfile(const std::filesystem::path &p) noexcept(false)
    -> bool {
  log<LogEvent::message>("Trying to load the library file {:s}", p.string());

  // Obtain file status and determine it's properties.
  std::error_code ec;
  const auto status = std::filesystem::status(p, ec);

  if (ec.value() != 0) {
    log<LogEvent::error>("OS error while trying to open {:s}: {:s}", p.string(), ec.message());
    throw std::runtime_error("IO error");
  }

  if (!std::filesystem::exists(status)) {
    log<LogEvent::error>("The file {:s} does not exist.", p.string());
    throw std::runtime_error("IO error");
  }

  if (!std::filesystem::is_regular_file(status)) {
    log<LogEvent::error>("The file {:s} is not regular.", p.string());
    throw std::runtime_error("IO error");
  }

  if (p.extension().string() != ".so") {
    log<LogEvent::error>("The file {:s} is not a .so file", p.string());
    throw std::runtime_error("IO error");
  }

  return true;
}

auto grlensing::shared_library::load_lib(const std::filesystem::path &p) noexcept(false)
    -> libhandle_t {
  auto lib = libhandle_t(dlopen(p.c_str(), RTLD_NOW), dlclose);

  if (lib == nullptr) {
    log<LogEvent::error>("Error while loading shared library: {:s}", dlerror());
    throw std::runtime_error("IO error");
  }

  return lib;
}