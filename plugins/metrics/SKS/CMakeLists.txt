
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# -----------------------------------------
# 1) Project
# -----------------------------------------

project(
  grlensing_SKS_metric
  VERSION 1.0.0
  LANGUAGES CXX
)

# -----------------------------------------
# 2) Sources
# -----------------------------------------

# Option file
configure_file(
  "${grlensing_SKS_metric_SOURCE_DIR}/include/options_in.txt"
  "${grlensing_SKS_metric_SOURCE_DIR}/include/options.hpp"
)

# Core source files
set(
  grlensing_SKS_metric_HEADER_LIST
  "${grlensing_SKS_metric_SOURCE_DIR}/include/api_macros.hpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/include/aux_functions.hpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/include/options.hpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/include/SKS.hpp"
  "${GRLensing_SOURCE_DIR}/include/kernel.hpp"
  "${GRLensing_SOURCE_DIR}/include/tensor_types.hpp"
)

  set(
  grlensing_SKS_metric_SOURCE_LIST
  "${grlensing_SKS_metric_SOURCE_DIR}/source/sx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/sy.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dsx_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dsy_dt.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/d2sx_dt2.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/d2sy_dt2.cpp"
  
  "${grlensing_SKS_metric_SOURCE_DIR}/source/T.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/X.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/Y.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/Z.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dT_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dT_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dT_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dT_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dX_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dX_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dX_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dX_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dY_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dY_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dY_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dY_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dZ_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dZ_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dZ_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dZ_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/g00.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g01.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g02.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g03.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g11.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g12.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g13.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g22.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g23.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/g33.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg00_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg00_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg00_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg00_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg01_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg01_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg01_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg01_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg02_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg02_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg02_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg02_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg03_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg03_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg03_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg03_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg11_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg11_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg11_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg11_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg12_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg12_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg12_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg12_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg13_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg13_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg13_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg13_dz.cpp"

  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg22_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg22_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg22_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg22_dz.cpp"
  
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg23_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg23_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg23_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg23_dz.cpp"
  
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg33_dt.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg33_dx.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg33_dy.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/dg33_dz.cpp"
  
  "${grlensing_SKS_metric_SOURCE_DIR}/source/grad_lapse.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/grad_ushift.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/l_shift.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/lapse.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/ll_extrinsic.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/ll_smetric.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/SKS.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/spatial_christoffel.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/u_shift.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/ul_extrinsic.cpp"
  "${grlensing_SKS_metric_SOURCE_DIR}/source/uu_smetric.cpp"
)

# Plugin module library
add_library(grlensing_SKS_metric MODULE ${grlensing_SKS_metric_HEADER_LIST} ${grlensing_SKS_metric_SOURCE_LIST})
target_compile_features(grlensing_SKS_metric PRIVATE cxx_std_20)
set_target_properties(grlensing_SKS_metric PROPERTIES OUTPUT_NAME "grlensing_SKS_metric")

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
    target_compile_options(grlensing_SKS_metric PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(grlensing_SKS_metric PUBLIC -g3 -O0 -Wall -Wextra -Werror -pedantic -pedantic-errors -fno-omit-frame-pointer -fsanitize=address,undefined)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_SKS_metric PUBLIC /Wall /Wextra /permissive- /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Generating a release build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_SKS_metric PUBLIC -O2 -Wall -Wextra -pedantic)
    target_link_options(grlensing_SKS_metric PUBLIC -O2 -Wall -Wextra -pedantic)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_SKS_metric PUBLIC /O2)
  endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Profile")
  message(STATUS "Generating a profiling build system")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(grlensing_SKS_metric PUBLIC -g3 -O2 -pg)
    target_link_options(grlensing_SKS_metric PUBLIC -g3 -O2 -pg)
  elseif(MSVC)
    # TODO: get the actual flags and dll libraryfor msvc
    target_compile_options(grlensing_SKS_metric PUBLIC /O2)
  endif()
endif()

# -----------------------------------------
# 4) Link dependencies
# -----------------------------------------

target_include_directories(
  grlensing_SKS_metric PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_link_libraries(grlensing_SKS_metric PRIVATE 
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
  CONAN_PKG::range-v3
  CONAN_PKG::fmt
  CONAN_PKG::docopt.cpp
  CONAN_PKG::yaml-cpp
  CONAN_PKG::sundials
  ${GRLENSING_DLL_LIBRARY}
)
