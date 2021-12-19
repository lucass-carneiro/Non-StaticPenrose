
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_isotropic_Schwarzschild_metric
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
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/include/options_in.txt"
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/include/options.hpp"
)

# Core source files
set(
  grlensing_isotropic_Schwarzschild_metric_HEADER_LIST
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/include/api_macros.hpp"
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/include/isotropic_Schwarzschild.hpp"
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/include/options.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
)

set(
  grlensing_isotropic_Schwarzschild_metric_SOURCE_LIST
  "${grlensing_isotropic_Schwarzschild_metric_SOURCE_DIR}/source/isotropic_Schwarzschild.cpp"
)

# Plugin module library
add_library(grlensing_isotropic_Schwarzschild_metric MODULE ${grlensing_isotropic_Schwarzschild_metric_HEADER_LIST} ${grlensing_isotropic_Schwarzschild_metric_SOURCE_LIST})
target_compile_features(grlensing_isotropic_Schwarzschild_metric PRIVATE cxx_std_20)
set_target_properties(grlensing_isotropic_Schwarzschild_metric PROPERTIES OUTPUT_NAME "grlensing_isotropic_Schwarzschild_metric")

# -----------------------------------------
# 3) Compilers flags and options
# -----------------------------------------

# Enforcing standard compliance and optimization flags on different compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(grlensing_isotropic_Schwarzschild_metric PUBLIC -Wall -Wextra -pedantic -O2)
elseif(MSVC)
  target_compile_options(grlensing_isotropic_Schwarzschild_metric PUBLIC /Wall /Wextra /permissive- /O2)
endif()

# -----------------------------------------
# 4) Link dependencies
# -----------------------------------------

target_include_directories(
  grlensing_isotropic_Schwarzschild_metric PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_link_libraries(grlensing_isotropic_Schwarzschild_metric PRIVATE yaml-cpp)
