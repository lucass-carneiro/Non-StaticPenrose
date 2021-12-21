#ifndef GRLENSING_TIME_INTEGRATION_HPP
#define GRLENSING_TIME_INTEGRATION_HPP

// clang-format off
#include "log.hpp"
#include "storage_server.hpp"
#include "metric_server.hpp"

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <tuple>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <yaml-cpp/node/node.h>
// clang-format on

namespace grlensing {

/**
 * Alias for a pointer to a file write robject
 */
using writer_ptr = std::unique_ptr<grlensing::storage_server::trajectory_writer>;

/**
 * Holds the necessary data to plot a single trajectory in a given spacetime.
 */
struct trajectory_config {
  realtype initial_time;
  realtype final_time;

  realtype initial_V1;
  realtype initial_V2;
  realtype initial_V3;

  realtype initial_X1;
  realtype initial_X2;
  realtype initial_X3;

  realtype initial_EN;

  realtype background_radius;
  realtype energy_threshold;

  int output_times;

  int particle_type;

  /**
   * Parses a trajectory configuration from a yaml file.
   *
   * @param trajectory_file The yaml file contain data for the trajectory.
   * @return A new trajectory_config object.
   */
  trajectory_config(const YAML::Node &trajectory_file)
      : initial_time(trajectory_file["initial_time"].as<realtype>()),
        final_time(trajectory_file["final_time"].as<realtype>()),
        initial_V1(trajectory_file["initial_V1"].as<realtype>()),
        initial_V2(trajectory_file["initial_V2"].as<realtype>()),
        initial_V3(trajectory_file["initial_V3"].as<realtype>()),
        initial_X1(trajectory_file["initial_X1"].as<realtype>()),
        initial_X2(trajectory_file["initial_X2"].as<realtype>()),
        initial_X3(trajectory_file["initial_X3"].as<realtype>()),
        initial_EN(trajectory_file["initial_EN"].as<realtype>()),
        background_radius(trajectory_file["background_radius"].as<realtype>()),
        energy_threshold(trajectory_file["energy_threshold"].as<realtype>()),
        output_times(trajectory_file["output_times"].as<int>()),
        particle_type(trajectory_file["particle_type"].as<int>()) {}
};

/**
 * The type of the user data necessary during integration
 */
using user_data_t
    = std::tuple<const metric_server::metric_ptr &, const grlensing::trajectory_config &>;

/**
 * Contains general settings of the integrator
 */
struct integrator_config {
  realtype absolute_tolerance;
  realtype relative_tolerance;
  realtype convergence_coefficient;

  int max_error_test_fails;
  int max_nonlinear_iterations;
  int max_steps;
  int predictor_method_order;

  /**
   * Construct integrator settings from configuration file.
   *
   * @param config_file The yaml config file with the integrator configuration.
   * @return A new integrator config object.
   */
  integrator_config(const YAML::Node &config_file)
      : absolute_tolerance(config_file["integrator_settings"]["absolute_tolerance"].as<realtype>()),
        relative_tolerance(config_file["integrator_settings"]["relative_tolerance"].as<realtype>()),
        convergence_coefficient(
            config_file["integrator_settings"]["convergence_coefficient"].as<realtype>()),
        max_error_test_fails(config_file["integrator_settings"]["max_error_test_fails"].as<int>()),
        max_nonlinear_iterations(
            config_file["integrator_settings"]["max_nonlinear_iterations"].as<int>()),
        max_steps(config_file["integrator_settings"]["max_steps"].as<int>()),
        predictor_method_order(
            config_file["integrator_settings"]["predictor_method_order"].as<int>()) {

    if (max_error_test_fails < 0) {
      log<LogEvent::warning>("The parameter \"max_error_test_fails\" is negative. It will be "
                             "automatically turned into a positive value.");
      max_error_test_fails = -1 * max_error_test_fails;
    }

    if (max_nonlinear_iterations < 0) {
      log<LogEvent::warning>("The parameter \"max_nonlinear_iterations\" is negative. It will be "
                             "automatically turned into a positive value.");
      max_nonlinear_iterations = -1 * max_nonlinear_iterations;
    }

    if (max_steps < 0) {
      log<LogEvent::warning>("The parameter \"max_steps\" is negative. It will be "
                             "automatically turned into a positive value.");
      max_steps = -1 * max_steps;
    }

    if (predictor_method_order < 0) {
      log<LogEvent::warning>("The parameter \"predictor_method_order\" is negative. It will be "
                             "automatically turned into a positive value.");
      predictor_method_order = -1 * predictor_method_order;
    }
  }
};

/**
 * The geodesic equation system in ADM form.
 *
 * @param t The time parameter of the ODE system. In this case, the coordinate time.
 * @param y The left hand side (state vector) of the system.
 * @param ydot The right hand side of the ODE system.
 * @param user_data Pointer to user data necessary to be used in the integration.
 * @return A status code resulting from the process of writing the rhs.
 */
auto adm_geodesic_system(realtype t, N_Vector y, N_Vector ydot, void *user_data) noexcept -> int;

/**
 * The jacobian of the geodesic equation system in ADM form.
 *
 * @param t The time parameter of the ODE system. In this case, the coordinate time.
 * @param y The left hand side (state vector) of the system.
 * @param fy TODO: What is this parameter?
 * @param J The Jacobian matrix to write.
 * @param user_data Pointer to user data necessary to be used in the integration.
 * @param tmp1 TODO: What is this?
 * @param tmp2 TODO: What is this?
 * @param tmp3 TODO: What is this?
 * @return A status code resulting from the process of writing the Jacobian.
 */
auto adm_geodesic_jacobian(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept -> int;
/**
 * Detects if a particle has colided with the background during the integration
 * of the system or if it was swalloed by a black hole.
 *
 * Background collision occurs when a particle trajectory intersects with a sphere of user supplied
 * radius. Swallowing by the astrophysical object occures when the energy becomes greater then a
 * user supplied threshold.
 *
 * @param t The time parameter of the ODE system. In this case, the coordinate time.
 * @param y The left hand side (state vector) of the system.
 * @param cs A pointer for a vector containing the LHS of the equation system that determines if a
 * colision ocurred.
 * @param user_data Pointer to user data necessary to be used in the integration.
 * @return A status code resulting from the process of writing the Jacobian.
 */
auto detect_collisions(realtype t, N_Vector y, realtype *cs, void *user_data) noexcept -> int;

/**
 * Integrates the geodesic equation.
 *
 * @param int_conf The integrator configuration data.
 * @param writer A file writer that will be used for outputting the data.
 * @param metric A spacetime background metric.
 * @param traj_conf The trajectory configuration data.
 */
void integrate(const integrator_config &int_conf, const writer_ptr &writer,
               const metric_server::metric_ptr &metric,
               const trajectory_config &traj_conf) noexcept(false);

/**
 * Check the return values of SUNDIAL routines for errors.
 *
 * @param flagvalue The value of  the flag to check.
 * @param funcname The name of the function being checked.
 * @param opt Wha to to actually check:
 * 0 Checks for returned null pointers in SUNDIALS functions.
 * 1 Checks the returned integer values from SUNDIALS functions
 */
template <typename T> void check_flag(const T &flagvalue, const char *funcname) noexcept(false) {
  if constexpr (std::is_pointer<T>::value) {
    if (flagvalue == nullptr) {
      log<LogEvent::error>("SUNDIALS error. The function {:s} returned a null pointer.", funcname);
      throw std::runtime_error("SUNDIALS null pointer error.");
    }
  } else if constexpr (std::is_integral<T>::value) {
    if (flagvalue < 0) {
      log<LogEvent::error>("SUNDIALS error. The function {:s} returned {:d}.", funcname, flagvalue);
      throw std::runtime_error("SUNDIALS function failure.");
    }
  }
}

/**
 * Normalizes velocities.
 *
 * When the user chooses the initial velocitis for a particle, these might not be physicall, that
 * is, their norm might be greater than one. If the user chooses to compute the trajectory of
 * photons, however, the velocity norm must be exaclty one. This function normalizes the user
 * supplied velocities in this case to ensure that the propper normalization is achieved.
 *
 * @param config A trajectory configuration data.
 * @param metric A spacetime metric objject.
 */
void normalize(trajectory_config &config, const metric_server::metric_ptr &metric);

} // namespace grlensing

#endif // GRLENSING_TIME_INTEGRATION_HPP
