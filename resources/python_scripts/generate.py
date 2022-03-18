# generate.py
# Generates boilerplate code for a new metric

import os
import subprocess

# The name of the metric
name = "isotropic_Schwarzschild"
display_name = "Isotropic Schwarzschild"

# If the metric is static, the time coordinate will not appear explicitly
# on the ADM functions, which avoid compiler warnings and errors.
is_time_dependand = False
is_space_dependant = False

# The spacetime coordinates
coords = ["t", "x", "y", "z"]

#-------------------------------------------------------------------------------
# Configuration header template
#-------------------------------------------------------------------------------

options_template = f"""
#ifndef GRLENSING_{name.upper()}_METRIC_OPTIONS_HPP
#define GRLENSING_{name.upper()}_METRIC_OPTIONS_HPP

#define GRLENSING_{name.upper()}_METRIC_VERSION_MAJOR @grlensing_{name}_metric_VERSION_MAJOR@ // NOLINT
#define GRLENSING_{name.upper()}_METRIC_VERSION_MINOR @grlensing_{name}_metric_VERSION_MINOR@ // NOLINT
#define GRLENSING_{name.upper()}_METRIC_VERSION_PATCH @grlensing_{name}_metric_VERSION_PATCH@ // NOLINT

#define GRLENSING_{name.upper()}_METRIC_SYSTEM_@CMAKE_SYSTEM_NAME@
#define GRLENSING_{name.upper()}_METRIC_COMPILER_@CMAKE_CXX_COMPILER_ID@

#endif // GRLENSING_{name.upper()}_METRIC_OPTIONS_HPP
"""

#-------------------------------------------------------------------------------
# API macros header template
#-------------------------------------------------------------------------------

api_macros_template = f"""
#ifndef GRLENSING_{name.upper()}_METRIC_API_MACROS_HPP
#define GRLENSING_{name.upper()}_METRIC_API_MACROS_HPP

#include "../../../../include/kernel.hpp"
#include "options.hpp"

#ifdef GRLENSING_{name.upper()}_METRIC_COMPILER_MSVC
#  if defined(GRLENSING_{name.upper()}_METRIC_STATICLIB)
#    define GRLENSING_{name.upper()}_METRIC_API
#  else
#    ifdef GRLENSING_{name.upper()}_METRIC_SOURCE
#      define GRLENSING_{name.upper()}_METRIC_API __declspec(dllexport)
#    else
#      define GRLENSING_{name.upper()}_METRIC_API __declspec(dllimport)
#    endif
#  endif

#elif defined(GRLENSING_{name.upper()}_METRIC_COMPILER_GNU)
#  if defined(GRLENSING_{name.upper()}_METRIC_STATICLIB)
#    define GRLENSING_{name.upper()}_METRIC_API
#  else
#    if defined(GRLENSING_{name.upper()}_METRIC_SOURCE)
#      define GRLENSING_{name.upper()}_METRIC_API __attribute__((visibility("default")))
#    else
// NOLINTNEXTLINE
#      define GRLENSING_{name.upper()}_METRIC_API __attribute__((visibility("default")))
#    endif
#  endif

#else

#  error Unknown compiler, please implement shared library API macros

#endif

#endif // GRLENSING_{name.upper()}_METRIC_API_MACROS_HPP
"""

#-------------------------------------------------------------------------------
# CMakeLists.txt
#-------------------------------------------------------------------------------

cmakelists = f"""
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_{name}_metric
  VERSION 1.0.0
  LANGUAGES CXX
)

# -----------------------------------------
# 2) Dependencies
# -----------------------------------------

# YAML parser
find_package(yaml-cpp REQUIRED)

# -----------------------------------------
# 3) Sources
# -----------------------------------------

# Option file
configure_file(
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/include/options_in.txt"
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/include/options.hpp"
)

# Core source files
set(
  grlensing_{name}_metric_HEADER_LIST
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/include/api_macros.hpp"
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/include/{name}.hpp"
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/include/options.hpp"
  "$\u007bGRLensing_SOURCE_DIR\u007d/include/kernel.hpp"
)

set(
  grlensing_{name}_metric_SOURCE_LIST
  "$\u007bgrlensing_{name}_metric_SOURCE_DIR\u007d/source/{name}.cpp"
)

# Plugin module library
add_library(grlensing_{name}_metric MODULE $\u007bgrlensing_{name}_metric_HEADER_LIST\u007d $\u007bgrlensing_{name}_metric_SOURCE_LIST\u007d)
target_compile_features(grlensing_{name}_metric PRIVATE cxx_std_20)
set_target_properties(grlensing_{name}_metric PROPERTIES OUTPUT_NAME "grlensing_{name}_metric")

# -----------------------------------------
# 3) Compilers flags and options
# -----------------------------------------

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(GRLENSING_DLL_LIBRARY dl)
elseif(MSVC)
  # TODO: get the actual flags and dll libraryfor msvc
  set(GRLENSING_DLL_LIBRARY TODO.dll)
endif()

if($<CONFIG:Debug>)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(GRLensing PUBLIC -g -Wall -Wextra -Werror -pedantic -pedantic-errors -O2 -fsanitize=address,undefined)
    target_link_options(GRLensing PUBLIC -g -Wall -Wextra -Werror -pedantic -pedantic-errors -O2 -fsanitize=address,undefined)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(GRLensing PUBLIC /Wall /Wextra /permissive- /O2)
  endif()
elseif($<CONFIG:Release>)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(GRLensing PUBLIC -O2)
    target_link_options(GRLensing PUBLIC -O2)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(GRLensing PUBLIC /O2)
  endif()
endif()

# -----------------------------------------
# 4) Link dependencies
# -----------------------------------------

target_include_directories(
  grlensing_{name}_metric PUBLIC $<BUILD_INTERFACE:$\u007bPROJECT_SOURCE_DIR\u007d/include>
  $<INSTALL_INTERFACE:include/$\u007bPROJECT_NAME\u007d-$\u007bPROJECT_VERSION\u007d>
)

target_link_libraries(grlensing_{name}_metric PRIVATE yaml-cpp)
"""

#-------------------------------------------------------------------------------
# Main header file
#-------------------------------------------------------------------------------

main_header = f"""
#ifndef GRLENSING_{name.upper()}_METRIC_HPP
#define GRLENSING_{name.upper()}_METRIC_HPP

#include "api_macros.hpp"

namespace grlensing \u007b

class {name} final : public metric_server::adm_metric \u007b
public:
  GRLENSING_{name.upper()}_METRIC_API void load_parameters(const YAML::Node &) final;

  GRLENSING_{name.upper()}_METRIC_API auto lapse(double, double, double, double)
      -> double final;

  GRLENSING_{name.upper()}_METRIC_API auto l_shift(double, double, double, double)
      -> metric_server::spatial_vector final;

  GRLENSING_{name.upper()}_METRIC_API auto u_shift(double, double, double, double) -> metric_server::spatial_vector final;

  GRLENSING_{name.upper()}_METRIC_API auto ll_smetric(double, double, double, double) -> metric_server::spatial_matrix final;

  GRLENSING_{name.upper()}_METRIC_API auto uu_smetric(double, double, double, double) -> metric_server::spatial_matrix final;

  GRLENSING_{name.upper()}_METRIC_API auto ll_extrinsic(double, double, double, double) -> metric_server::spatial_matrix final;

  GRLENSING_{name.upper()}_METRIC_API auto ul_extrinsic(double, double, double, double) -> metric_server::spatial_matrix final;

  GRLENSING_{name.upper()}_METRIC_API auto spatial_christoffel(double, double, double, double) -> metric_server::chirstofell_t final;

  GRLENSING_{name.upper()}_METRIC_API auto grad_lapse(double, double, double, double) -> metric_server::spatial_vector final;

  GRLENSING_{name.upper()}_METRIC_API auto grad_ushift(double, double, double, double) -> metric_server::spatial_matrix final;

  GRLENSING_{name.upper()}_METRIC_API auto name() -> std::string_view final;
\u007d;

extern "C" GRLENSING_{name.upper()}_METRIC_API auto get_engine_version() -> unsigned;

extern "C" GRLENSING_{name.upper()}_METRIC_API void register_plugin(kernel &kernel);

\u007d // namespace grlensing

#endif // GRLENSING_{name.upper()}_METRIC_HPP ";
"""

#-------------------------------------------------------------------------------
# Main source file
#-------------------------------------------------------------------------------

if(is_time_dependand and is_space_dependant):
    function_signature = f"double {coords[0]}, double {coords[1]}, double {coords[2]}, double {coords[3]}"
elif(is_time_dependand and not is_space_dependant):
    function_signature = f"double {coords[0]}, double, double, double"
elif(not is_time_dependand and is_space_dependant):
    function_signature = f"double, double {coords[1]}, double {coords[2]}, double {coords[3]}"
else:
    function_signature = f"double, double, double, double"


main_source = f"""
#include \"{name}.hpp\"

using grlensing::kernel;
using grlensing::metric_server;
using grlensing::{name};

GRLENSING_{name.upper()}_METRIC_API void {name}::load_parameters(const YAML::Node &) \u007b \u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::lapse({function_signature}) -> double \u007b
  return 1.0;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::l_shift({function_signature}) -> metric_server::spatial_vector \u007b
  return metric_server::spatial_vector \u007b 0.0, 0.0, 0.0 \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::u_shift({function_signature}) -> metric_server::spatial_vector \u007b
  return metric_server::spatial_vector \u007b 0.0, 0.0, 0.0 \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::ll_smetric({function_signature}) -> metric_server::spatial_matrix \u007b
  return metric_server::spatial_matrix \u007b \u007b \u007b 1.0, 0.0, 0.0 \u007d, \u007b 0.0, 1.0, 0.0 \u007d, \u007b 0.0, 0.0, 1.0 \u007d \u007d \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::uu_smetric({function_signature}) -> metric_server::spatial_matrix \u007b
  return metric_server::spatial_matrix \u007b \u007b \u007b 1.0, 0.0, 0.0 \u007d, \u007b 0.0, 1.0, 0.0 \u007d, \u007b 0.0, 0.0, 1.0 \u007d \u007d \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::ll_extrinsic({function_signature}) -> metric_server::spatial_matrix \u007b
  return metric_server::spatial_matrix \u007b \u007b \u007b 0.0, 0.0, 0.0 \u007d, \u007b  0.0, 0.0, 0\u007d, \u007b  0.0, 0.0, 0.0 \u007d \u007d \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::ul_extrinsic({function_signature}) -> metric_server::spatial_matrix \u007b
  return metric_server::spatial_matrix \u007b \u007b \u007b 0.0, 0.0, 0.0 \u007d, \u007b  0.0, 0.0, 0\u007d, \u007b  0.0, 0.0, 0.0 \u007d \u007d \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::spatial_christoffel({function_signature}) -> metric_server::chirstofell_t \u007b
  metric_server::chirstofell_t Gamma;
  
  Gamma[0][0][0] = 0.0;
  Gamma[0][0][1] = 0.0;
  Gamma[0][0][2] = 0.0;
  Gamma[0][1][1] = 0.0;
  Gamma[0][1][2] = 0.0;
  Gamma[0][2][2] = 0.0;

  Gamma[1][0][0] = 0.0;
  Gamma[1][0][1] = 0.0;
  Gamma[1][0][2] = 0.0;
  Gamma[1][1][1] = 0.0;
  Gamma[1][1][2] = 0.0;
  Gamma[1][2][2] = 0.0;

  Gamma[2][0][0] = 0.0;
  Gamma[2][0][1] = 0.0;
  Gamma[2][0][2] = 0.0;
  Gamma[2][1][1] = 0.0;
  Gamma[2][1][2] = 0.0;
  Gamma[2][2][2] = 0.0;

  Gamma[0][1][0] = Gamma[0][0][1];
  Gamma[0][2][0] = Gamma[0][0][2];
  Gamma[0][2][1] = Gamma[0][1][2];
  Gamma[1][1][0] = Gamma[1][0][1];
  Gamma[1][2][0] = Gamma[1][0][2];
  Gamma[1][2][1] = Gamma[1][1][2];
  Gamma[2][1][0] = Gamma[2][0][1];
  Gamma[2][2][0] = Gamma[2][0][2];
  Gamma[2][2][1] = Gamma[2][1][2];

  return Gamma;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::grad_lapse({function_signature}) -> metric_server::spatial_vector \u007b
  return metric_server::spatial_vector \u007b 0.0, 0.0, 0.0 \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::grad_ushift(double, double, double, double) -> metric_server::spatial_matrix \u007b
  return metric_server::spatial_matrix \u007b \u007b \u007b 0.0, 0.0, 0.0 \u007d, \u007b  0.0, 0.0, 0\u007d, \u007b  0.0, 0.0, 0.0 \u007d \u007d \u007d;
\u007d

GRLENSING_{name.upper()}_METRIC_API auto {name}::name() -> std::string_view \u007b
  return "{display_name}";
\u007d

extern "C" GRLENSING_{name.upper()}_METRIC_API auto get_engine_version() -> unsigned \u007b
  return unsigned(1);
\u007d

extern "C" GRLENSING_{name.upper()}_METRIC_API void register_plugin(kernel &kernel) \u007b
  kernel.get_metric_server().add_metric(
      std::unique_ptr<metric_server::adm_metric>(new {name}()));
\u007d

"""

#-------------------------------------------------------------------------------
# Generation
#-------------------------------------------------------------------------------

os.makedirs(name + "/source", exist_ok=True)
os.makedirs(name + "/include", exist_ok=True)

file = open(name + "/include/api_macros.hpp", mode="w")
file.write(api_macros_template)
file.close()

file = open(name + f"/include/{name}.hpp", mode="w")
file.write(main_header)
file.close()

file = open(name + f"/include/options_in.txt", mode="w")
file.write(options_template)
file.close()

file = open(name + f"/source/{name}.cpp", mode="w")
file.write(main_source)
file.close()

file = open(name + "/CMakeLists.txt", mode="w")
file.write(cmakelists)
file.close()

subprocess.run(
  [
    "clang-format",
    "-i",
    name + f"/include/{name}.hpp",
    name + f"/source/{name}.cpp"
  ]
)