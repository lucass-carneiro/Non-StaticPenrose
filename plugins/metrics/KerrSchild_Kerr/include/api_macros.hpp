
#ifndef GRLENSING_KERRSCHILD_KERR_METRIC_API_MACROS_HPP
#define GRLENSING_KERRSCHILD_KERR_METRIC_API_MACROS_HPP

#include "../../../../include/kernel.hpp"
#include "options.hpp"

#ifdef GRLENSING_KERRSCHILD_KERR_METRIC_COMPILER_MSVC
#  if defined(GRLENSING_KERRSCHILD_KERR_METRIC_STATICLIB)
#    define GRLENSING_KERRSCHILD_KERR_METRIC_API
#  else
#    ifdef GRLENSING_KERRSCHILD_KERR_METRIC_SOURCE
#      define GRLENSING_KERRSCHILD_KERR_METRIC_API __declspec(dllexport)
#    else
#      define GRLENSING_KERRSCHILD_KERR_METRIC_API __declspec(dllimport)
#    endif
#  endif

#elif defined(GRLENSING_KERRSCHILD_KERR_METRIC_COMPILER_GNU)
#  if defined(GRLENSING_KERRSCHILD_KERR_METRIC_STATICLIB)
#    define GRLENSING_KERRSCHILD_KERR_METRIC_API
#  else
#    if defined(GRLENSING_KERRSCHILD_KERR_METRIC_SOURCE)
#      define GRLENSING_KERRSCHILD_KERR_METRIC_API __attribute__((visibility("default")))
#    else
// NOLINTNEXTLINE
#      define GRLENSING_KERRSCHILD_KERR_METRIC_API __attribute__((visibility("default")))
#    endif
#  endif

#else

#  error Unknown compiler, please implement shared library API macros

#endif

#endif // GRLENSING_KERRSCHILD_KERR_METRIC_API_MACROS_HPP
