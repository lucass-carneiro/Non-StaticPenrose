
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_KerrSchild_Kerr_metric
  VERSION 1.0.0
  LANGUAGES CXX
)

# -----------------------------------------
# 2) Dependencies
# -----------------------------------------

# YAML parser
find_package(yaml-cpp REQUIRED)
find_package(fmt REQUIRED)

# -----------------------------------------
# 3) Sources
# -----------------------------------------

# Option file
configure_file(
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/options_in.txt"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/options.hpp"
)

# Core source files
set(
  grlensing_KerrSchild_Kerr_metric_HEADER_LIST
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/api_macros.hpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/aux_functions.hpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/KerrSchild_Kerr.hpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/include/options.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
)

set(
  grlensing_KerrSchild_Kerr_metric_SOURCE_LIST
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/aux_functions.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/christoffel.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/grad_lapse.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/grad_ushift.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/KerrSchild_Kerr.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/ll_extrinsic.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/metric.cpp"
  "${grlensing_KerrSchild_Kerr_metric_SOURCE_DIR}/source/ul_extrinsic.cpp"
)

# Plugin module library
add_library(grlensing_KerrSchild_Kerr_metric MODULE ${grlensing_KerrSchild_Kerr_metric_HEADER_LIST} ${grlensing_KerrSchild_Kerr_metric_SOURCE_LIST})
target_compile_features(grlensing_KerrSchild_Kerr_metric PRIVATE cxx_std_20)
set_target_properties(grlensing_KerrSchild_Kerr_metric PROPERTIES OUTPUT_NAME "grlensing_KerrSchild_Kerr_metric")

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
  grlensing_KerrSchild_Kerr_metric PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_link_libraries(grlensing_KerrSchild_Kerr_metric PRIVATE yaml-cpp fmt)
