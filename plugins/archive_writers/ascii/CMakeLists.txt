cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_std_ascii_writer_plugin
  VERSION 1.0.0
  LANGUAGES CXX
)

# -----------------------------------------
# 2) Sources
# -----------------------------------------

# Option file
configure_file(
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/include/options_in.txt"
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/include/options.hpp"
)

# Core source files
set(
  grlensing_std_ascii_writer_plugin_HEADER_LIST
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/include/api_macros.hpp"
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/include/options.hpp"
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/include/ascii_writer.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
)

set(
  grlensing_std_ascii_writer_plugin_SOURCE_LIST
  "${grlensing_std_ascii_writer_plugin_SOURCE_DIR}/source/ascii_writer.cpp"
)

# Plugin module library
add_library(grlensing_std_ascii_writer_plugin MODULE 
    ${grlensing_std_ascii_writer_plugin_HEADER_LIST}
    ${grlensing_std_ascii_writer_plugin_SOURCE_LIST}
)
target_compile_features(grlensing_std_ascii_writer_plugin PRIVATE cxx_std_20)
set_target_properties(grlensing_std_ascii_writer_plugin PROPERTIES OUTPUT_NAME "grlensing_ascii_writer")

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
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(grlensing_std_ascii_writer_plugin PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC /Wall /Wextra /permissive- /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Generating a release build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC -O2 -Wall -Wextra -pedantic)
    target_link_options(grlensing_std_ascii_writer_plugin PUBLIC -O2 -Wall -Wextra -pedantic)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
  message(STATUS "Generating a profiling build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC -g3 -O2 -pg)
    target_link_options(grlensing_std_ascii_writer_plugin PUBLIC -g3 -O2 -pg)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_std_ascii_writer_plugin PUBLIC /O2)
  endif()
endif()

# -----------------------------------------
# 7) Link dependencies
# -----------------------------------------

target_link_libraries(grlensing_std_ascii_writer_plugin PRIVATE CONAN_PKG::fmt)

target_include_directories(
  grlensing_std_ascii_writer_plugin PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_link_libraries(grlensing_std_ascii_writer_plugin PRIVATE 
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
  CONAN_PKG::range-v3
  CONAN_PKG::fmt
  CONAN_PKG::docopt.cpp
  CONAN_PKG::yaml-cpp
  CONAN_PKG::sundials
  ${GRLENSING_DLL_LIBRARY}
)
