cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_std_zip_plugin
  VERSION 1.0.0
  LANGUAGES CXX
)

# -----------------------------------------
# 2) Sources
# -----------------------------------------

# Option file
configure_file(
  "${grlensing_std_zip_plugin_SOURCE_DIR}/include/options_in.txt"
  "${grlensing_std_zip_plugin_SOURCE_DIR}/include/options.hpp"
)

# Core source files
set(
  GRLENSING_STD_ZIP_PLUGIN_HEADER_LIST
  "${grlensing_std_zip_plugin_SOURCE_DIR}/include/api_macros.hpp"
  "${grlensing_std_zip_plugin_SOURCE_DIR}/include/options.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
)

set(
  GRLENSING_STD_ZIP_PLUGIN_SOURCE_LIST
  "${grlensing_std_zip_plugin_SOURCE_DIR}/source/zip_reader.cpp"
)

# Plugin module library
add_library(grlensing_std_zip_plugin MODULE ${GRLENSING_STD_ZIP_PLUGIN_HEADER_LIST} ${GRLENSING_STD_ZIP_PLUGIN_SOURCE_LIST})
target_compile_features(grlensing_std_zip_plugin PRIVATE cxx_std_20)
set_target_properties(grlensing_std_zip_plugin PROPERTIES OUTPUT_NAME "grlensing_zip_archive")

# -----------------------------------------
# 3) Compilers flags and options
# -----------------------------------------

# Enforcing standard compliance and optimization flags on different compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(grlensing_std_zip_plugin PUBLIC -Wall -Wextra -pedantic -O2)
elseif(MSVC)
  target_compile_options(grlensing_std_zip_plugin PUBLIC /Wall /Wextra /permissive- /O2)
endif()

# -----------------------------------------
# 7) Link dependencies
# -----------------------------------------

target_include_directories(
  grlensing_std_zip_plugin PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)