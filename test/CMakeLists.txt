set(GRLENSING_TESTS_SOURCE_LIST ${GRLENSING_SOURCE_LIST})
list(REMOVE_ITEM GRLENSING_TESTS_SOURCE_LIST ${GRLensing_SOURCE_DIR}/source/main.cpp)
list(APPEND GRLENSING_TESTS_SOURCE_LIST "${GRLensing_SOURCE_DIR}/test/tests.cpp")

enable_testing()

add_executable(tests ${GRLENSING_HEADER_LIST} ${GRLENSING_TESTS_SOURCE_LIST})
target_compile_features(tests PRIVATE cxx_std_20)

target_include_directories(
    tests PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

target_link_libraries(tests PRIVATE 
  MPI::MPI_CXX
  OpenMP::OpenMP_CXX
  CONAN_PKG::range-v3
  CONAN_PKG::fmt
  CONAN_PKG::yaml-cpp
  CONAN_PKG::sundials
  CONAN_PKG::gtest 
  ${GRLENSING_DLL_LIBRARY}
)

include(GoogleTest)
gtest_discover_tests(tests)