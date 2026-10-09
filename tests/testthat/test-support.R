
test_that("Control presets", {
  ## control() is a lowercase alias for the Control() constructor
  expect_inherits(control(), "Control")
  expect_equal(control(), Control())

  ## Defaults are the pragmatic, fast-ish settings used for most runs
  expect_equal(Control()$ode_tol_rel, 1e-4)
  expect_equal(Control()$ode_tol_abs, 1e-4)
  expect_equal(Control()$ode_step_size_max, 5)
  expect_equal(Control()$node_gradient_direction, -1)
  expect_equal(Control()$schedule_eps, 2e-2)

  ## control_accurate() tightens the ODE and schedule tolerances
  acc <- control_accurate()
  expect_inherits(acc, "Control")
  expect_equal(acc$ode_tol_rel, 1e-6)
  expect_equal(acc$ode_tol_abs, 1e-6)
  expect_equal(acc$ode_step_size_max, 1e-1)
  expect_equal(acc$schedule_eps, 1e-3)
})

test_that("control_tf24() sets TF24's step control and keeps the rest", {
  tf <- control_tf24()
  expect_inherits(tf, "Control")
  expect_equal(tf$ode_tol_rel, 3e-5)
  expect_equal(tf$ode_tol_abs, 3e-9)
  expect_equal(tf$ode_tol_factor_soil, 10)
  expect_equal(tf$ode_tol_factor_max, 100)
  expect_equal(tf$ode_step_size_max, 15 / 365)
  expect_equal(tf$ode_soil_substep_max_uptake, 0.1)
  ## The absolute tolerance follows the relative one.
  expect_equal(control_tf24(1e-5)$ode_tol_abs, 1e-9)

  ## Every other field is the base's.
  base <- Control(ode_split_sign_changes = TRUE, schedule_eps = 1e-3)
  set <- c("ode_tol_rel", "ode_tol_abs", "ode_tol_factor_soil",
           "ode_tol_factor_max", "ode_step_size_max", "ode_soil_substep_max_uptake")
  rest <- setdiff(names(base), set)
  expect_identical(unclass(control_tf24(base = base))[rest],
                   unclass(base)[rest])
})
