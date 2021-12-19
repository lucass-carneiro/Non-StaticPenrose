#ifndef GRLENSING_API_MACROS_HPP
#define GRLENSING_API_MACROS_HPP

#include "options.hpp"

#ifdef GRLENSING_COMPILER_MSVC
#  if defined(GRLENSING_STATICLIB)
#    define GRLENSING_API
#  else
#    ifdef GRLENSING_SOURCE
#      define GRLENSING_API __declspec(dllexport)
#    else
#      define GRLENSING_API __declspec(dllimport)
#    endif
#  endif

#elif defined(GRLENSING_COMPILER_GNU)

#  if defined(GRLENSING_STATICLIB)
#    define GRLENSING_API
#  else
#    if defined(GRLENSING_SOURCE)
#      define GRLENSING_API __attribute__((visibility("default")))
#    else
#      define GRLENSING_API __attribute__((visibility("default"))) // NOLINT
#    endif
#  endif

#else

#  error Unknown compiler, please implement shared library API macros

#endif

#endif // GRLENSING_API_MACROS_HPP