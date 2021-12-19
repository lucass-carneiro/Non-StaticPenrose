#ifndef GRLENSING_ASCII_WRITER_API_MACROS_HPP
#define GRLENSING_ASCII_WRITER_API_MACROS_HPP

#include "../../../../include/kernel.hpp"
#include "options.hpp"

#ifdef GRLENSING_ASCII_WRITER_COMPILER_MSVC
#  if defined(GRLENSING_ASCII_WRITER_STATICLIB)
#    define GRLENSING_ASCII_WRITER_API
#  else
#    ifdef GRLENSING_ASCII_WRITER_SOURCE
#      define GRLENSING_ASCII_WRITER_API __declspec(dllexport)
#    else
#      define GRLENSING_ASCII_WRITER_API __declspec(dllimport)
#    endif
#  endif

#elif defined(GRLENSING_ASCII_WRITER_COMPILER_GNU)
#  if defined(GRLENSING_ASCII_WRITER_STATICLIB)
#    define GRLENSING_ASCII_WRITER_API
#  else
#    if defined(GRLENSING_ASCII_WRITER_SOURCE)
#      define GRLENSING_ASCII_WRITER_API __attribute__((visibility("default")))
#    else
#      define GRLENSING_ASCII_WRITER_API __attribute__((visibility("default")))
#    endif
#  endif

#else

#  error Unknown compiler, please implement shared library API macros

#endif

#endif // GRLENSING_ASCII_WRITER_API_MACROS_HPP
