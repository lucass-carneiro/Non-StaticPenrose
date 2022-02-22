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

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(GRLENSING_DLL_LIBRARY dl)
elseif(MSVC)
  # TODO: get the actual flags and dll libraryfor msvc
  set(GRLENSING_DLL_LIBRARY TODO.dll)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Generating a debug build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_std_zip_plugin PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(grlensing_std_zip_plugin PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_zip_plugin PUBLIC /Wall /Wextra /permissive- /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Generating a release build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_std_zip_plugin PUBLIC -O2 -Wall -Wextra -pedantic)
    target_link_options(grlensing_std_zip_plugin PUBLIC -O2 -Wall -Wextra -pedantic)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_zip_plugin PUBLIC /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
  message(STATUS "Generating a profiling build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_std_zip_plugin PUBLIC -g3 -O2 -pg)
    target_link_options(grlensing_std_zip_plugin PUBLIC -g3 -O2 -pg)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_zip_plugin PUBLIC /O2)
  endif()
endif()

# -----------------------------------------
# 7) Link dependencies
# -----------------------------------------

target_include_directories(
  grlensing_std_zip_plugin PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)