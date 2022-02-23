set(GRLENSING_TESTS_SOURCE_LIST
  "${GRLensing_SOURCE_DIR}/test/tests.cpp"
)

enable_testing()

add_executable(tests ${GRLENSING_TESTS_SOURCE_LIST})
target_compile_features(tests PRIVATE cxx_std_20)

target_link_libraries(tests PRIVATE CONAN_PKG::gtest MPI::MPI_CXX)
include(GoogleTest)
gtest_discover_tests(tests)