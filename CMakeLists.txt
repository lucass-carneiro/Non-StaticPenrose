cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  GRLensing
  VERSION 1.0.0
  LANGUAGES CXX
  )

# Generate compile_commands.json for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

# -----------------------------------------
# 2) In-source build guard
# -----------------------------------------

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
    "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()

# -----------------------------------------
# 3) External dependencies
# -----------------------------------------

# The CMake package manager
include(cmake/CPM.cmake)

# Point to the external folder for packages
list(APPEND CMAKE_MODULE_PATH 
  "${CMAKE_CURRENT_LIST_DIR}/external/lib/cmake/"
  "${CMAKE_CURRENT_LIST_DIR}/external/share/"
)

set(docopt_DIR "${CMAKE_CURRENT_LIST_DIR}/external/lib/cmake/docopt")
set(SUNDIALS_DIR "${CMAKE_CURRENT_LIST_DIR}/external/lib/cmake/sundials")

# Message passing for distributed parallelism
find_package(MPI REQUIRED)

# Command line arguemnt parsing
find_package(docopt REQUIRED)

# YAML parser
find_package(yaml-cpp REQUIRED)

# Printf like formating
find_package(fmt REQUIRED)

# SUNDIALS ODE integrators
find_package(SUNDIALS REQUIRED)

# -----------------------------------------
# 4) Source files
# -----------------------------------------

# Available options
option(GRLENSING_USE_LOG_COLOR "Use colors on log outputs" ON)

# Option file
configure_file("${GRLensing_SOURCE_DIR}/include/options_in.txt" "${GRLensing_SOURCE_DIR}/include/options.hpp")

# Core source files
set(
  GRLENSING_HEADER_LIST
  "${GRLensing_SOURCE_DIR}/include/api_macros.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
  "${GRLensing_SOURCE_DIR}/include/log.hpp"
  "${GRLensing_SOURCE_DIR}/include/metric_server.hpp"
  "${GRLensing_SOURCE_DIR}/include/options.hpp"
  "${GRLensing_SOURCE_DIR}/include/plugin.hpp"
  "${GRLensing_SOURCE_DIR}/include/shared_library.hpp"
  "${GRLensing_SOURCE_DIR}/include/storage_server.hpp"
  "${GRLensing_SOURCE_DIR}/include/time_integration.hpp"
)

set(
  GRLENSING_SOURCE_LIST
  "${GRLensing_SOURCE_DIR}/source/kernel.cpp"
  "${GRLensing_SOURCE_DIR}/source/main.cpp"
  "${GRLensing_SOURCE_DIR}/source/metric_server.cpp"
  "${GRLensing_SOURCE_DIR}/source/plugin.cpp"
  "${GRLensing_SOURCE_DIR}/source/shared_library.cpp"
  "${GRLensing_SOURCE_DIR}/source/storage_server.cpp"
  "${GRLensing_SOURCE_DIR}/source/adm_geodesic_callbacks.cpp"
  "${GRLensing_SOURCE_DIR}/source/adm_geodesic_equations.cpp"
  "${GRLensing_SOURCE_DIR}/source/adm_geodesic_integrate.cpp"
  "${GRLensing_SOURCE_DIR}/source/adm_geodesic_jacobian.cpp"
)

# Executable target
add_executable(GRLensing ${GRLENSING_HEADER_LIST} ${GRLENSING_SOURCE_LIST})
target_compile_features(GRLensing PRIVATE cxx_std_20)
set_target_properties(GRLensing PROPERTIES OUTPUT_NAME "grlensing")

# -----------------------------------------
# 5) Compilers flags and options
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
# 7) Link dependencies
# -----------------------------------------

target_link_libraries(GRLensing PRIVATE 
  MPI::MPI_CXX 
  ${GRLensing_SOURCE_DIR}/external/lib/libfmt.a #For some reason, CMake can't find fmt
  docopt 
  yaml-cpp 
  SUNDIALS::arkode 
  ${GRLENSING_DLL_LIBRARY}
)

target_include_directories(
  GRLensing PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

# -----------------------------------------
# 8) Compile modules
# -----------------------------------------

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/plugins/archive_readers)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/plugins/metrics)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/plugins/archive_writers)
