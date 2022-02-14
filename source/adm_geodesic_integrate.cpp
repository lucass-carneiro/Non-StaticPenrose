// clang-format off
#include "log.hpp"
#include "metric_server.hpp"
#include "time_integration.hpp"

#include <cstddef>
#include <fmt/printf.h>

#include <memory>
#include <cstdio>

#include <arkode/arkode_arkstep.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sunlinsol/sunlinsol_dense.h>
// clang-format on

// https://github.com/LLNL/sundials/blob/master/examples/arkode/C_serial/ark_robertson_root.c
// https://raw.githubusercontent.com/LLNL/sundials/master/doc/arkode/ark_examples.pdf
// https://raw.githubusercontent.com/LLNL/sundials/master/doc/arkode/ark_guide.pdf
void grlensing::integrate(const integrator_config &int_conf, const writer_ptr &writer,
                          const metric_server::metric_ptr &metric,
                          const trajectory_config &traj_conf) noexcept(false) {

  constexpr const auto system_size = 7;
  constexpr const auto callback_size = 2;

  const auto output_times = traj_conf.output_times;
  const auto output_time_step
      = std::abs(traj_conf.final_time - traj_conf.initial_time) / output_times;

  using nvec = _generic_N_Vector;
  using nvec_deleter = void (*)(N_Vector);

  using sunmat = _generic_SUNMatrix;
  using sunmat_deleter = void (*)(SUNMatrix);

  using sunlinsol = _generic_SUNLinearSolver;
  using sunlinsol_deleter = int (*)(SUNLinearSolver);

  using arkstep = void;
  using arksetp_deleter = void (*)(arkstep *);

  int flag = 0, rtflag = 0;
  std::array<int, callback_size> rootsfound = {0, 0};

  // Data structure initialization
  auto y = std::unique_ptr<nvec, nvec_deleter>(N_VNew_Serial(system_size), N_VDestroy);
  check_flag(y, "N_VNew_Serial");

  auto atols = std::unique_ptr<nvec, nvec_deleter>(N_VNew_Serial(system_size), N_VDestroy);
  check_flag(atols, "N_VNew_Serial");

  // Loading initial data
  NV_Ith_S(y.get(), 0) = traj_conf.initial_V1; // NOLINT
  NV_Ith_S(y.get(), 1) = traj_conf.initial_V2; // NOLINT
  NV_Ith_S(y.get(), 2) = traj_conf.initial_V3; // NOLINT
  NV_Ith_S(y.get(), 3) = traj_conf.initial_X1; // NOLINT
  NV_Ith_S(y.get(), 4) = traj_conf.initial_X2; // NOLINT
  NV_Ith_S(y.get(), 5) = traj_conf.initial_X3; // NOLINT
  NV_Ith_S(y.get(), 6) = traj_conf.initial_EN; // NOLINT

  /* Call ARKStepCreate to initialize the ARK timestepper module and
   * specify the right-hand side function in y'=f(t,y), the inital time
   * T0, and the initial dependent variable vector y.  Note: since this
   * problem is fully implicit, we set f_E to null and f_I to the system.
   */
  auto arkode_mem = std::unique_ptr<arkstep, arksetp_deleter>(
      ARKStepCreate(nullptr, adm_geodesic_system, traj_conf.initial_time, y.get()),
      [](arkstep *p) -> void { ARKStepFree(&p); });
  check_flag(arkode_mem, "ARKStepCreate");

  // Set tolerances
  for (std::size_t i = 0; i < system_size; i++)
    NV_Ith_S(atols, i) = int_conf.absolute_tolerance; // NOLINT

  // Increase max error test fails
  flag = ARKStepSetMaxErrTestFails(arkode_mem.get(), int_conf.max_error_test_fails);
  check_flag(flag, "ARKStepSetMaxErrTestFails");

  // Increase max nonlinear iterations
  flag = ARKStepSetMaxNonlinIters(arkode_mem.get(), int_conf.max_nonlinear_iterations);
  check_flag(flag, "ARKStepSetMaxNonlinIters");

  // Update nonlinear solver convergence coeff.
  flag = ARKStepSetNonlinConvCoef(arkode_mem.get(), int_conf.convergence_coefficient);
  check_flag(flag, "ARKStepSetNonlinConvCoef");

  // Increase max number of steps
  flag = ARKStepSetMaxNumSteps(arkode_mem.get(), int_conf.max_steps);
  check_flag(flag, "ARKStepSetMaxNumSteps");

  // Specify maximum-order predictor
  flag = ARKStepSetPredictorMethod(arkode_mem.get(), int_conf.predictor_method_order);
  check_flag(flag, "ARKStepSetPredictorMethod");

  // Specify tolerances
  flag = ARKStepSVtolerances(arkode_mem.get(), int_conf.relative_tolerance, atols.get());
  check_flag(flag, "ARKStepSStolerances");

  // Specify the root-finding function, having 2 equations
  flag = ARKStepRootInit(arkode_mem.get(), callback_size, detect_collisions);
  check_flag(flag, "ARKStepRootInit");

  // Initialize dense matrix data structure and solver
  auto A = std::unique_ptr<sunmat, sunmat_deleter>(SUNDenseMatrix(system_size, system_size),
                                                   SUNMatDestroy);
  check_flag(A, "SUNDenseMatrix");

  auto LS = std::unique_ptr<sunlinsol, sunlinsol_deleter>(SUNLinSol_Dense(y.get(), A.get()),
                                                          SUNLinSolFree);
  check_flag(A, "SUNLinSol_Dense");

  // Attach matrix and linear solver
  flag = ARKStepSetLinearSolver(arkode_mem.get(), LS.get(), A.get());
  check_flag(flag, "ARKStepSetLinearSolver");

  // Set the Jacobian routine
  flag = ARKStepSetJacFn(arkode_mem.get(), adm_geodesic_jacobian);
  check_flag(flag, "ARKStepSetJacFn");

  // Set user data
  user_data_t user_data = std::make_tuple(std::ref(metric), std::ref(traj_conf));
  flag = ARKStepSetUserData(arkode_mem.get(), static_cast<void *>(&user_data));
  check_flag(flag, "ARKStepSetUserData");

  // output initial condition to disk
  writer->push_real(traj_conf.initial_time);
  for (std::size_t i = 0; i < system_size; i++)
    writer->push_real(NV_Ith_S(y.get(), i)); // NOLINT

  metric_server::spatial_vector Vi
      = {NV_Ith_S(y.get(), 0), NV_Ith_S(y.get(), 1), NV_Ith_S(y.get(), 2)}; // NOLINT

  metric_server::spatial_vector Xi
      = {NV_Ith_S(y.get(), 3), NV_Ith_S(y.get(), 4), NV_Ith_S(y.get(), 5)}; // NOLINT

  // NOLINTNEXTLINE
  double El = NV_Ith_S(y.get(), 6);
  auto Eg = compute_global_energy(metric, traj_conf.initial_time, Vi, Xi, El);
  writer->push_final_real(Eg);

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
   * prints results.  Stops when the final time has been reached
   */
  auto t = traj_conf.initial_time;

  auto tout = t + output_time_step;
  auto iout = 0;

  while (true) {
    // call integrator
    flag = ARKStepEvolve(arkode_mem.get(), tout, y.get(), &t, ARK_NORMAL);
    check_flag(flag, "ARKStepEvolve");

    // access/print solution
    writer->push_real(t);
    for (std::size_t i = 0; i < system_size; i++)
      writer->push_real(NV_Ith_S(y.get(), i)); // NOLINT

    // NOLINTNEXTLINE
    Vi = {NV_Ith_S(y.get(), 0), NV_Ith_S(y.get(), 1), NV_Ith_S(y.get(), 2)};

    // NOLINTNEXTLINE
    Xi = {NV_Ith_S(y.get(), 3), NV_Ith_S(y.get(), 4), NV_Ith_S(y.get(), 5)};

    // NOLINTNEXTLINE
    El = NV_Ith_S(y.get(), 6);
    Eg = compute_global_energy(metric, t, Vi, Xi, El);
    writer->push_final_real(Eg);

    // Exit due to root found.
    if (flag == ARK_ROOT_RETURN) {
      rtflag = ARKStepGetRootInfo(arkode_mem.get(), rootsfound.data());
      check_flag(rtflag, "ARKStepGetRootInfo");
      writer->push_metadata("background colision: ",
                            rootsfound[0] ? std::size_t(1) : std::size_t(0));
      writer->push_metadata("swallowed: ", rootsfound[1] ? std::size_t(1) : std::size_t(0));
      break;
    }

    // Exit due to solver error.
    if (flag >= 0) {
      iout++;
      tout = t + output_time_step;
    } else {
      log<LogEvent::error>("Solver failure, stopping integration\n");
      break;
    }

    // Exit due to enough outpus performed.
    if (iout == output_times)
      break;
  }

  // Print some final statistics
  long int nst = 0, nst_a = 0, nfe = 0, nfi = 0, nsetups = 0;
  long int nje = 0, nfeLS = 0, nni = 0, ncfn = 0, netf = 0, nge = 0;

  flag = ARKStepGetNumSteps(arkode_mem.get(), &nst);
  check_flag(flag, "ARKStepGetNumSteps");

  flag = ARKStepGetNumStepAttempts(arkode_mem.get(), &nst_a);
  check_flag(flag, "ARKStepGetNumStepAttempts");

  flag = ARKStepGetNumRhsEvals(arkode_mem.get(), &nfe, &nfi);
  check_flag(flag, "ARKStepGetNumRhsEvals");

  flag = ARKStepGetNumLinSolvSetups(arkode_mem.get(), &nsetups);
  check_flag(flag, "ARKStepGetNumLinSolvSetups");

  flag = ARKStepGetNumErrTestFails(arkode_mem.get(), &netf);
  check_flag(flag, "ARKStepGetNumErrTestFails");

  flag = ARKStepGetNumNonlinSolvIters(arkode_mem.get(), &nni);
  check_flag(flag, "ARKStepGetNumNonlinSolvIters");

  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem.get(), &ncfn);
  check_flag(flag, "ARKStepGetNumNonlinSolvConvFails");

  flag = ARKStepGetNumJacEvals(arkode_mem.get(), &nje);
  check_flag(flag, "ARKStepGetNumJacEvals");

  flag = ARKStepGetNumLinRhsEvals(arkode_mem.get(), &nfeLS);
  check_flag(flag, "ARKStepGetNumLinRhsEvals");

  flag = ARKStepGetNumGEvals(arkode_mem.get(), &nge);
  check_flag(flag, "ARKStepGetNumGEvals");

  writer->push_metadata("   Internal solver steps", std::size_t(nst));
  writer->push_metadata("   Attempted Internal solver steps", std::size_t(nst_a));

  writer->push_metadata("   Total RHS evals Fe", std::size_t(nfe));
  writer->push_metadata("   Total RHS evals Fi", std::size_t(nfi));

  writer->push_metadata("   Total linear solver setups", std::size_t(nsetups));
  writer->push_metadata("   Total RHS evals for setting up the linear system", std::size_t(nfeLS));
  writer->push_metadata("   Total number of Jacobian evaluations", std::size_t(nje));
  writer->push_metadata("   Total number of Newton iterations", std::size_t(nni));
  writer->push_metadata("   Total root-function g evals", std::size_t(nge));
  writer->push_metadata("   Total number of nonlinear solver convergence failures",
                        std::size_t(ncfn));
  writer->push_metadata("   Total number of error test failures", std::size_t(netf));

  if (traj_conf.particle_type == 0) {
    writer->push_metadata("   V1 after normalization", traj_conf.initial_V1);
    writer->push_metadata("   V2 after normalization", traj_conf.initial_V2);
    writer->push_metadata("   V3 after normalization", traj_conf.initial_V3);
  }
}

void grlensing::normalize(trajectory_config &traj_conf, const metric_server::metric_ptr &metric) {
  // The current velocity norm
  const metric_server::spatial_vector V
      = {traj_conf.initial_V1, traj_conf.initial_V2, traj_conf.initial_V3};
  const auto gamma
      = metric->ll_smetric(0.0, traj_conf.initial_X1, traj_conf.initial_X2, traj_conf.initial_X3);

  const auto norm_V = gamma[0][0] * V[0] * V[0] + 2.0 * gamma[0][1] * V[0] * V[1]
                      + 2.0 * gamma[0][2] * V[0] * V[2] + gamma[1][1] * V[1] * V[1]
                      + 2.0 * gamma[1][2] * V[1] * V[2] + gamma[2][2] * V[2] * V[2];

  // For photons, $p_\mu p^\mu = 0$ which implies $V_i V^i = 1$
  if (traj_conf.particle_type == 0) {
    traj_conf.initial_V1 /= norm_V;
    traj_conf.initial_V2 /= norm_V;
    traj_conf.initial_V3 /= norm_V;
  } else if (!(norm_V < 1.0)) {
    log<LogEvent::error>(
        "The velocity values used represent a massive particle that move with superluminal speeds. "
        "Ensure that $gamma_{{ij}} V^i V^j < 1$. The current value is {0:.16e}",
        norm_V);
    throw std::runtime_error("Unphysical particle error");
  }
}