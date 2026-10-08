
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
  expect_equal(tf$ode_weight_soil, 10)
  expect_equal(tf$ode_weight_max, 100)
  expect_equal(tf$ode_step_size_max, 15 / 365)
  expect_equal(tf$ode_soil_alone_share, 0.1)
  ## The absolute tolerance follows the relative one.
  expect_equal(control_tf24(1e-5)$ode_tol_abs, 1e-9)

  ## Every other field is the base's.
  base <- Control(ode_split_sign_changes = TRUE, schedule_eps = 1e-3)
  set <- c("ode_tol_rel", "ode_tol_abs", "ode_weight_soil", "ode_weight_max",
           "ode_step_size_max", "ode_soil_alone_share")
  rest <- setdiff(names(base), set)
  expect_identical(unclass(control_tf24(base = base))[rest],
                   unclass(base)[rest])
})

test_that("control_window() weighs each step by the offspring still to be earned", {
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
  earned <- vapply(pilot$history, function(h) {
    o <- h$species[[1]]$net_reproduction_ratio_by_node
    sum(weight[seq_along(o)] * o)
  }, 0)
  R <- 1 - earned / earned[length(earned)]
  at <- vapply(pilot$history, function(h) h$time, 0)
  factor_at <- function(w, t) w$ode_weight_factors[findInterval(t, w$ode_weight_times)]

  win <- control_window(pilot, base = ctrl)
  expect_equal(win$ode_weight_times[1], 0)
  expect_false(is.unsorted(win$ode_weight_times))
  expect_equal(factor_at(win, at), 1 / pmin(pmax(R / 0.1, 0.01), 1))
  expect_true(any(factor_at(win, at) > 1 & factor_at(win, at) < 100))
  expect_equal(range(win$ode_weight_factors), c(1, 100))
  expect_equal(max(control_window(pilot, base = ctrl, r_min = 0.1)$ode_weight_factors), 10)
  rest <- setdiff(names(ctrl), c("ode_weight_times", "ode_weight_factors"))
  expect_identical(unclass(win)[rest], unclass(ctrl)[rest])

  ## An invader that matures later earns later, so the window that protects it
  ## weighs no step more and some less.
  guarded <- control_window(pilot, list(stand(3)), base = ctrl)
  expect_true(all(guarded$ode_weight_factors <= win$ode_weight_factors))
  expect_true(any(guarded$ode_weight_factors < win$ode_weight_factors))

  height <- run_scm(stand(2), env(), control_tf24(1e-3))
  expect_error(control_window(height), "birth-date coordinate")
})

test_that("diagnose_scm() estimates a run's node error from coarser runs", {
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 8
  times <- seq(0, 6.4, length.out = 9)
  stand <- function(hmat, lma = 0.1978791, at = times) {
    p <- add_strategies(p0, trait_matrix(c(lma, hmat), c("lma", "hmat")))
    p$node_schedule_times <- list(at)
    p
  }
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 2)
  ctrl <- control_tf24(1e-3, Control(node_density_in_birth_date = TRUE))
  pulse <- function(p) events(events_default(p), rainfall_pulse(time = 1.5, depth = 0))
  ## lma at a million gives a newborn shorter than one leaf segment, so its walk
  ## throws.
  d <- diagnose_scm(stand(2), env, ctrl, list(late = stand(3), broken = stand(2, 1e6)),
                    events = pulse(stand(2)))
  expect_equal(d$nodes, c(9, 5, 3))

  ## The coarser run is the run at every other introduction, its pulse kept.
  ln_J <- function(at) {
    scm <- run_scm(stand(2, at = at), env, ctrl, events = pulse(stand(2, at = at)))
    log(sum(scm$offspring_production))
  }
  q <- d$quantities
  resident <- q[q$run == "resident" & q$quantity == "ln J", ]
  expect_identical(resident$half, ln_J(times[c(1, 3, 5, 7, 9)]))
  ## Its correction moves the run toward one at twice the introductions.
  finer <- ln_J(seq(0, 6.4, length.out = 17))
  expect_lt(abs(resident$value + resident$error - finer), abs(resident$value - finer))
  expect_true(all(c("lma", "hmat") %in% q$quantity[q$run == "late"]))
  expect_true(all(is.finite(q$value[q$run %in% c("resident", "late")])))

  late <- d$distance[d$distance$run == "late", ]
  expect_equal(late$parameter, "hmat")
  expect_equal(late$distance, log(1.5))

  ## A walk that throws is recorded at each run, and the rest still answer.
  expect_false("broken" %in% q$run)
  expect_equal(d$failures$run, rep("broken", 3))
  expect_equal(d$failures$nodes, c(9, 5, 3))
  expect_equal(unique(d$failures$stage), "walk")
  expect_match(d$failures$message[1], "L_tip")

  alone <- diagnose_scm(stand(2), env, ctrl, gradient = FALSE)
  expect_equal(alone$quantities$quantity, "ln J")
  expect_equal(nrow(alone$failures), 0)

  expect_error(diagnose_scm(stand(2), env, ctrl, list(stand(3))), "names each invader")
  expect_error(diagnose_scm(stand(2, at = times[1:3]), env, ctrl), "four introductions")
  two <- add_strategies(stand(2), trait_matrix(0.3, "lma"))
  expect_error(diagnose_scm(two, env, ctrl), "one species")
})
