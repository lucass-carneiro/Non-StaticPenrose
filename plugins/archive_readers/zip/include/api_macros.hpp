#ifndef GRLENSING_ZIP_PLUGIN_CONFIG_HPP
#define GRLENSING_ZIP_PLUGIN_CONFIG_HPP

#include "../../../../include/kernel.hpp"
#include "options.hpp"

#ifdef GRLENSING_ZIP_PLUGIN_COMPILER_MSVC
#  if defined(GRLENSING_ZIP_PLUGIN_STATICLIB)
#    define GRLENSING_ZIP_PLUGIN_API
#  else
#    ifdef GRLENSING_ZIP_PLUGIN_SOURCE
#      define GRLENSING_ZIP_PLUGIN_API __declspec(dllexport)
#    else
#      define GRLENSING_ZIP_PLUGIN_API __declspec(dllimport)
#    endif
#  endif

#elif defined(GRLENSING_ZIP_PLUGIN_COMPILER_GNU)

#  if defined(GRLENSING_ZIP_PLUGIN_STATICLIB)
#    define GRLENSING_ZIP_PLUGIN_API
#  else
#    if defined(GRLENSING_ZIP_PLUGIN_SOURCE)
#      define GRLENSING_ZIP_PLUGIN_API __attribute__((visibility("default")))
#    else
#      define GRLENSING_ZIP_PLUGIN_API __attribute__((visibility("default")))
#    endif
#  endif

#else

#  error Unknown compiler, please implement shared library API macros

#endif

#endif // GRLENSING_ZIP_PLUGIN_CONFIG_HPP
