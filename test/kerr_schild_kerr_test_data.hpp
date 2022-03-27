#ifndef GRLENSING_TESTS_KERR_SCHILD_KERR_TEST_DATA_HPP
#define GRLENSING_TESTS_KERR_SCHILD_KERR_TEST_DATA_HPP

#include "tensor_types.hpp"

namespace grlensing_tests {

/**
 * The radius function in Kerr-Schild coordiantes of the Kerr metric.
 *
 * @param a The black hole dimensionlens spin.
 * @param x The x coordinate value.
 * @param y The y coordinate value.
 * @param z The z coordinate value.
 */
auto r(double a, double x, double y, double z) -> double;

/**
 * The H function in Kerr-Schild coordiantes of the Kerr metric.
 *
 * @param M The mass of the black hole.
 * @param a The black hole dimensionlens spin.
 * @param x The x coordinate value.
 * @param y The y coordinate value.
 * @param z The z coordinate value.
 */
auto H(double M, double a, double x, double y, double z) -> double;

/**
 * The l vector in Kerr-Schild coordiantes of the Kerr metric.
 *
 * @param a The black hole dimensionlens spin.
 * @param x The x coordinate value.
 * @param y The y coordinate value.
 * @param z The z coordinate value.
 */
auto l(double a, double x, double y, double z) -> grlensing::callable_array<4, double>;

/**
 * The covariant 4-metric tensor in Kerr-Schild coordiantes of the Kerr metric.
 *
 * @param M The mass of the black hole.
 * @param a The black hole dimensionlens spin.
 * @param x The x coordinate value.
 * @param y The y coordinate value.
 * @param z The z coordinate value.
 */
auto ll_g(double M, double a, double x, double y, double z)
    -> grlensing::callable_matrix<4, 4, double>;

/**
 * The contravariant 4-metric tensor in Kerr-Schild coordiantes of the Kerr metric.
 *
 * @param M The mass of the black hole.
 * @param a The black hole dimensionlens spin.
 * @param x The x coordinate value.
 * @param y The y coordinate value.
 * @param z The z coordinate value.
 */
auto uu_g(double M, double a, double x, double y, double z)
    -> grlensing::callable_matrix<4, 4, double>;

} // namespace grlensing_tests

#endif // GRLENSING_TESTS_KERR_SCHILD_KERR_TEST_DATA_HPP