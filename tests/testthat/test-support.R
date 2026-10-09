
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

test_that("control_window() loosens the steps after most offspring is earned", {
  ## A stand that matures at 2 m, so that eight years hold its whole window.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 8
  stand <- function(hmat) {
    p <- add_strategies(p0, trait_matrix(c(0.1978791, hmat), c("lma", "hmat")))
    p$node_schedule_times <- list(seq(0, 6.4, length.out = 12))
    p
  }
  env <- function() {
    e <- Environment("TF24")
    e$extrinsic_drivers_set_constant("rainfall", 2)
    e
  }
  ctrl <- control_tf24(1e-3, Control(node_density_in_birth_date = TRUE))
  pilot <- SCM("TF24", "TF24_Env")(stand(2), env(), empty_events(), ctrl)
  pilot$collect <- TRUE
  pilot$record_trajectory <- TRUE
  pilot$run()

  ## The share still to be earned at the end of each interval, read off the
  ## nodes the run kept there, each weighted as offspring production weights it.
  sp <- pilot$patch$species[[1]]
  weight <- head(sp$establishment_weights, -1) * sp$patch_densities
  expect_equal(sum(weight * sp$net_reproduction_ratio_by_node) *
                 stand(2)$strategies[[1]]$pars$S_D,
               sum(pilot$offspring_production))
  ## The offspring produced by each step's end ends at the run's.
  produced <- pilot$offspring_produced_at_ode_times
  expect_length(produced, length(pilot$ode_times))
  expect_false(is.unsorted(produced))
  expect_equal(produced[length(produced)], sum(pilot$offspring_production))
  earned <- vapply(pilot$history, function(h) {
    o <- h$species[[1]]$net_reproduction_ratio_by_node
    sum(weight[seq_along(o)] * o)
  }, 0)
  left <- 1 - earned / earned[length(earned)]
  at <- vapply(pilot$history, function(h) h$time, 0)
  factor_at <- function(w, t) {
    w$ode_tol_factor_values[findInterval(t, w$ode_tol_factor_times)]
  }

  win <- control_window(pilot, base = ctrl)
  expect_equal(win$ode_tol_factor_times[1], 0)
  expect_false(is.unsorted(win$ode_tol_factor_times, strictly = TRUE))
  expect_equal(factor_at(win, at), pmin(pmax(0.1 / left, 1), 100))
  expect_true(any(factor_at(win, at) > 1 & factor_at(win, at) < 100))
  expect_equal(range(win$ode_tol_factor_values), c(1, 100))
  coarse <- control_window(pilot, base = ctrl, factor_limit = 10)
  expect_equal(max(coarse$ode_tol_factor_values), 10)
  rest <- setdiff(names(ctrl),
                  c("ode_tol_factor_times", "ode_tol_factor_values"))
  expect_identical(unclass(win)[rest], unclass(ctrl)[rest])

  ## An invader that matures later earns later, so the window that protects it
  ## loosens no step more and some less.
  guarded <- control_window(pilot, list(stand(3)), base = ctrl)
  expect_true(all(guarded$ode_tol_factor_values <= win$ode_tol_factor_values))
  expect_true(any(guarded$ode_tol_factor_values < win$ode_tol_factor_values))
  ## The pilot now holds the invader's walk, so it serves no second window.
  expect_true(pilot$invaded)
  expect_error(control_window(pilot, base = ctrl), "no invader has walked")

  height <- run_scm(stand(2), env(), control_tf24(1e-3))
  expect_error(control_window(height), "birth-date coordinate")
})
