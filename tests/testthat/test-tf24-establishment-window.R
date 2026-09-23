# TF24's newborns establish on their gate averaged over establishment_window
# rather than on the gate at the instant they are born. The average is one ODE
# state per species, after its nodes, and the boundary node is seated at it.

# One species per leaf mass per area given.
window_patch <- function(lma = 0.0825) {
  p <- scm_base_parameters("TF24")
  p <- add_strategies(p, trait_matrix(lma, "lma"),
                      birth_rate = rep(list(1), length(lma)))
  Patch("TF24", "TF24_Env")(p, Environment("TF24"),
                            Control(node_density_in_birth_date = TRUE))
}

# The gate a birth-size individual of `strategy` reads in `environment`.
instant_gate <- function(strategy, environment) {
  Individual("TF24", "TF24_Env")(strategy)$establishment_probability(environment)
}

test_that("the window is a TF24 parameter with a gradient column", {
  s <- TF24_Strategy()
  expect_equal(s$pars$establishment_window, 0.05)
  # It divides the rate the average relaxes at, so no window of zero or less
  # prepares.
  for (bad in c(0, -0.05, NaN, Inf)) {
    s$pars$establishment_window <- bad
    expect_error(Individual("TF24", "TF24_Env")(s), "establishment_window")
  }
  patch <- window_patch()
  expect_true("1.establishment_window" %in% ladder_trait_names_tf24(patch))
  expect_false("establishment_window" %in% names(census_undifferentiable_tf24()))
})

test_that("each species carries its averaged gate after its nodes", {
  patch <- window_patch(c(0.0825, 0.10))
  n_env <- patch$environment$ode_size
  node_size <- Node("TF24", "TF24_Env")(TF24_Strategy())$ode_size
  expect_equal(patch$ode_size, 2 + n_env)
  expect_equal(patch$node_ode_size, 0)

  patch$introduce_new_node(1, 0)
  patch$introduce_new_node(1, 0.01)
  patch$introduce_new_node(2, 0)
  expect_equal(patch$node_ode_size, 3 * node_size)
  expect_equal(patch$ode_size, 3 * node_size + 2 + n_env)
  sizes <- vapply(patch$species, function(s) s$ode_size, numeric(1))
  expect_equal(sizes, c(2, 1) * node_size + 1)

  # The flat state is each species' block in turn, then the environment's.
  y <- patch$ode_state
  expect_identical(y, c(patch$species[[1]]$ode_state, patch$species[[2]]$ode_state,
                        utils::tail(y, n_env)))
  gates <- c(y[[2 * node_size + 1]], y[[3 * node_size + 2]])
  expect_true(all(gates > 0 & gates <= 1))

  # FF16 establishes at the instant and carries nothing beside its nodes.
  pf <- scm_base_parameters("FF16")
  pf <- add_strategies(pf, trait_matrix(0.0825, "lma"))
  ff16 <- Patch("FF16", "FF16_Env")(pf, Environment("FF16"),
                                   Control(node_density_in_birth_date = TRUE))
  expect_equal(ff16$ode_size, 0)
  expect_equal(ff16$species[[1]]$ode_size, 0)
})

test_that("a patch that starts empty starts each window at its newborn's gate", {
  patch <- window_patch(c(0.0825, 0.10))
  p <- patch$parameters
  for (i in 1:2) {
    expect_equal(patch$species[[i]]$ode_state,
                 instant_gate(p$strategies[[i]], patch$environment))
  }
  # In equilibrium, so the rate is exactly zero.
  expect_identical(patch$ode_rates[1:2], c(0, 0))
})

test_that("held in one environment, the average relaxes to the gate on the window", {
  # With every other entry of the state held, the gate is one number and the
  # average obeys dE/dt = (gate - E)/window, whose solution from E0 is
  # gate + (E0 - gate) exp(-t/window).
  patch <- window_patch()
  window <- patch$parameters$strategies[[1]]$pars$establishment_window
  gate <- instant_gate(patch$parameters$strategies[[1]], patch$environment)
  y <- patch$ode_state
  rate <- function(e) {
    y[[1]] <- e
    patch$derivs(y, 0)[[1]]
  }
  for (e in c(0.2, 0.5, 1)) {
    expect_equal(rate(e), (gate - e) / window, tolerance = 1e-12)
  }

  # Integrated by the classical fourth-order scheme at a hundredth of the window.
  e0 <- 0.2
  h <- window / 100
  e <- e0
  for (k in seq_len(300)) {
    k1 <- rate(e)
    k2 <- rate(e + h / 2 * k1)
    k3 <- rate(e + h / 2 * k2)
    k4 <- rate(e + h * k3)
    e <- e + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  expect_equal(e, gate + (e0 - gate) * exp(-3), tolerance = 1e-9)
})

test_that("a node is seeded at the average, not at the gate of its instant", {
  # Light that swings deep enough to close the gate every year, so the average
  # and the gate part company between introductions.
  p <- ladder_parameters("fast", lifetime = 1.5)
  p$node_schedule_times <- list(c(0, 0.3, 0.55, 0.8, 1.05, 1.3))
  env <- Environment("TF24")
  t <- seq(0, 1.5, length.out = 400)
  env$extrinsic_drivers_set_variable("PPFD", t, 1800 * (1 + 0.97 * sin(2 * pi * t)))
  scm <- SCM("TF24", "TF24_Env")(p, env, empty_events(), ladder_control())
  scm$record_trajectory <- TRUE
  scm$run()
  trajectory <- scm$store_trajectory()
  names_node <- Node("TF24", "TF24_Env")(TF24_Strategy())$ode_names
  stride <- length(names_node)
  i_mortality <- match("mortality", names_node)
  i_density <- match("log_density", names_node)
  birth_rate <- p$strategies[[1]]$birth_rate_y

  seeded <- 0L
  for (r in trajectory[vapply(trajectory, function(r) r$introduction, TRUE)]) {
    n <- (length(r$state) - scm$patch$environment$ode_size - 1L) / stride
    newest <- (n - 1L) * stride
    average <- r$state[[n * stride + 1L]]
    expect_equal(r$state[[newest + i_mortality]], -log(average), tolerance = 1e-14)
    expect_equal(r$state[[newest + i_density]], log(birth_rate * average),
                 tolerance = 1e-14)
    seeded <- seeded + 1L
  }
  expect_equal(seeded, length(p$node_schedule_times[[1]]))

  # Non-vacuity: the average ranged over most of the gate's span.
  averages <- vapply(trajectory, function(r)
    r$state[[length(r$state) - scm$patch$environment$ode_size]], numeric(1))
  expect_lt(min(averages), 0.1)
  expect_gt(max(averages), 0.99)
})

test_that("the window's column agrees with a difference of whole runs", {
  # Both sides on one time grid, so the difference sees the model and not the
  # solver's choice of steps.
  build <- function(window) {
    p <- ladder_parameters("fast", lifetime = 1.5)
    p$node_schedule_times <- list(c(0, 0.3, 0.55, 0.8, 1.05, 1.3))
    p$strategies[[1]]$pars$establishment_window <- window
    env <- Environment("TF24")
    t <- seq(0, 1.5, length.out = 400)
    env$extrinsic_drivers_set_variable("PPFD", t,
                                       1800 * (1 + 0.97 * sin(2 * pi * t)))
    list(p = p, env = env)
  }
  window <- 0.05
  base <- build(window)
  scm <- SCM("TF24", "TF24_Env")(base$p, base$env, empty_events(), ladder_control())
  scm$run()
  result <- stand_gradient(scm)
  got <- result$gradient[, "1.establishment_window"]
  # The forward tangent of the same recording shares no transpose with the sweep.
  # Measured, they agree to 5e-16 of the column.
  tangent <- ladder_trajectory_tangent(
    scm, ladder_trait_direction(colnames(result$gradient),
                                "1.establishment_window"))$tangent
  expect_lt(max(abs(got - tangent)) / max(abs(got)), 1e-12)
  times <- scm$ode_times
  census_at <- function(w) {
    b <- build(w)
    b$p$ode_times <- times
    stand_census(ladder_run(b$p, env = b$env))
  }
  h <- window * 1e-4
  reference <- (census_at(window + h) - census_at(window - h)) / (2 * h)
  scale <- max(abs(reference))
  expect_gt(scale, 0)
  # Measured 4.8e-7, the difference's own floor at this step.
  expect_lt(max(abs(got - reference[names(got)])) / scale, 1e-5)
})
