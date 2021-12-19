#ifndef GRLENSING_SHARED_LIBRARY_HPP
#define GRLENSING_SHARED_LIBRARY_HPP

#include "api_macros.hpp"
#include "log.hpp"

#include <filesystem>
#include <memory>
#include <stdexcept>

#ifdef GRLENSING_SYSTEM_Linux
#  include <dlfcn.h>
#endif

namespace grlensing {

#ifdef GRLENSING_SYSTEM_Linux

/**
 * Represent and manage the lifetime of a shaderd (or dynamic) library
 */
class shared_library {
public:
  /**
   * Create a shared library object and load it.
   *
   * @param p Path to the dynamic library.
   */
  GRLENSING_API shared_library(const std::filesystem::path &p) noexcept(false)
      : handle((validate_libfile(p), load_lib(p))), name(p.filename().string()) {}

/* The use of reinterpret_cast from a void* (returned by dlsym) to a function pointer is
 * implementation defined behaviour. This is by design and unavoidable, however for systems that do
 * provide dlsym, this behaviour must be documented for conforming C++ compilers, thus we make sure
 * to be on a Linux system with a compiler that supports the conversion.
 */
#  if defined GRLENSING_SYSTEM_Linux                                                               \
      && (defined GRLENSING_COMPILER_GNU || defined GRLENSING_COMPILER_GNU)

  /**
   * Get a function pointer from a symbol inside the library.
   *
   * @param function_name The name the function to recover.
   * @return A function pointer to the requested symbol name.
   */
  template <typename T> auto get_function_pointer(const std::string &function_name) noexcept(false)
      -> T * {

    /* As per https://man7.org/linux/man-pages/man3/dlsym.3.html,
     * 1. Call dlerror() once to clear old error conditions.
     */
    dlerror();

    // 2. Attempt to get symbol with dlsym
    auto function_address = dlsym(handle.get(), function_name.c_str());

    // 3. call dlerror() again and save it's value.
    const auto *error = dlerror();

    // 4. Check the value from dlerror
    if (error != nullptr) {
      function_address = nullptr;
      log<LogEvent::error>("Function {:s} not found in {:s}", function_name, name);
      throw std::runtime_error("Dynamic library lookup error");
    } else {
      return reinterpret_cast<T *>(function_address); // NOLINT
    }
  }

#  else
#    error                                                                                         \
        "Dynamic function retrival from shared libraries not implemented for the current system/compiler."
#  endif

private:
  using libhandle_t = std::shared_ptr<void>;
  libhandle_t handle;
  const std::string name;

  /**
   * Validade a candidate library file.
   *
   * @param p Path  to a candidate library file.
   * @return true (false) for a valid (invalid) file
   */
  auto validate_libfile(const std::filesystem::path &p) noexcept(false) -> bool;

  /**
   * Load a shared library
   *
   * @param p Path to a library file.
   * @return Handle to the oppened library file
   */
  auto load_lib(const std::filesystem::path &p) noexcept(false) -> libhandle_t;
};

#endif

} // namespace grlensing

#endif // GRLENSING_SHARED_LIBRARY_HPP