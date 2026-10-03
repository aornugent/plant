# The SCM under ode_method = "ark": ARK4(3)6L[2]SA steps TF24's soil drainage and
# infiltration implicitly and everything else explicitly. A run replays on its
# own steps exactly and converges on Cash-Karp's answer; its sweep, for the stand
# and for an invader walking its field, matches central differences taken on its
# steps, through the soil as well as the traits.

# A two-node stand under seasonal rain, on the birth-date coordinate the sweep
# transposes, with the soil weight the long-drought runs use.
ark_fixture <- function() {
  p <- ladder_parameters("fast", lifetime = 1)
  p$node_schedule_times <- list(c(0, 0.37))
  list(p = p,
       env = ladder_environment(rain = 1, amplitude = 1, lifetime = 1),
       ctrl = ladder_control(ode_method = "ark", ode_weight_soil = 100))
}

# q with one parameter of its first strategy set past the hyperparameters, which
# is the parameter a gradient column differentiates.
ark_with_parameter <- function(q, name, value) {
  strategies <- q$strategies
  pars <- strategies[[1]]$pars
  pars[[name]] <- value
  strategies[[1]]$pars <- pars
  q$strategies <- strategies
  q
}

# The fixture's parameters pinned to a run's steps.
ark_on_steps <- function(p, run) {
  p$ode_times <- run$ode_times
  p$ode_step_sizes <- run$ode_step_sizes
  p
}

ark_offspring <- function(scm) stand_census(scm)[["offspring_production"]]

test_that("ARK converges on Cash-Karp's answer, stepping its own way", {
  f <- ark_fixture()
  ark <- run_scm(f$p, f$env, f$ctrl)
  ck <- run_scm(f$p, f$env, ladder_control(ode_weight_soil = 100))
  expect_false(identical(ark$ode_times, ck$ode_times))
  tight <- function(method) {
    run_scm(f$p, f$env, ladder_control(ode_method = method, ode_tol_rel = 1e-7,
                                       ode_tol_abs = 1e-7))
  }
  expect_equal(ark_offspring(tight("ark")), ark_offspring(tight("rkck")),
               tolerance = 1e-4)
})

test_that("an ARK run replayed on its own steps reproduces it", {
  # The walks and the sweep take a run's steps without its error estimate, so
  # the steps alone have to carry the run.
  f <- ark_fixture()
  run <- run_scm(f$p, f$env, f$ctrl)
  replayed <- run_scm(ark_on_steps(f$p, run), f$env, f$ctrl)
  expect_identical(replayed$ode_times, run$ode_times)
  expect_identical(replayed$ode_step_sizes, run$ode_step_sizes)
  expect_identical(replayed$patch$ode_state, run$patch$ode_state)
})

test_that("the ARK sweep matches central differences on the run's steps, through the traits and the soil", {
  f <- ark_fixture()
  stand <- run_scm(f$p, f$env, f$ctrl, record_trajectory = TRUE)
  traits <- c("lma", "a_dG2", "d_I")
  swept <- stand_gradient(stand, metrics = "offspring_production",
                          traits = paste0("1.", traits))$gradient[1, ]
  on_steps <- ark_on_steps(f$p, stand)
  differenced <- vapply(traits, function(name) {
    value <- f$p$strategies[[1]]$pars[[name]]
    h <- 1e-5 * value
    at <- function(x) {
      ark_offspring(run_scm(ark_with_parameter(on_steps, name, x), f$env, f$ctrl))
    }
    (at(value + h) - at(value - h)) / (2 * h)
  }, 0)
  # lma's difference has a floor near 1e-5; the other two agree to 2e-7.
  expect_equal(unname(swept / differenced), rep(1, 3), tolerance = 1e-4)

  # The soil, which every step carries through its implicit stages: the adjoint
  # the sweep ends holding at the first state, the soil alone there, against
  # replays from a perturbed soil. The top two layers carry what the stand reads
  # of it; the deeper ones sit below the difference's floor.
  metrics <- census_metric_names_tf24()
  at_first <- census_trait_gradient_tf24(stand)$at_first_state
  lambda <- at_first[[which(metrics == "offspring_production")]]
  base <- ladder_range_base_state(stand, 0L)
  expect_length(lambda, length(base))
  top <- length(base) - 9:8
  replayed <- function(state) {
    stats::setNames(ladder_census_initial_state_replay_tf24(stand, state, 0L),
                    metrics)[["offspring_production"]]
  }
  soil <- vapply(top, function(i) {
    h <- 1e-4 * base[i]
    (replayed(replace(base, i, base[i] + h)) -
     replayed(replace(base, i, base[i] - h))) / (2 * h)
  }, 0)
  expect_equal(lambda[top] / soil, rep(1, 2), tolerance = 1e-5)
})

test_that("an invader with the stand's traits walks an ARK run to its fitness, and its sweep matches differences of its walks", {
  # The walk takes the run's ARK steps in the fields the run recorded, so an
  # invader identical to the stand makes every evaluation the stand made.
  f <- ark_fixture()
  scm <- run_scm(f$p, f$env, f$ctrl)
  run_rr <- scm$net_reproduction_ratios
  expect_true(all(is.finite(run_rr)) && all(run_rr > 0))
  scm$run_mutant(f$p)
  expect_identical(scm$net_reproduction_ratios, run_rr)

  # Its selection gradient: the invasion kept no states, so the sweep repeats it
  # to keep them, and each side of a difference invades against the same
  # recording.
  traits <- c("lma", "a_dG2", "d_I")
  swept <- stand_gradient(scm, metrics = "offspring_production",
                          traits = paste0("1.", traits))$gradient[1, ]
  scm$record_trajectory <- FALSE
  differenced <- vapply(traits, function(name) {
    value <- f$p$strategies[[1]]$pars[[name]]
    h <- 1e-5 * value
    at <- function(x) {
      scm$run_mutant(ark_with_parameter(f$p, name, x))
      ark_offspring(scm)
    }
    (at(value + h) - at(value - h)) / (2 * h)
  }, 0)
  expect_equal(unname(swept / differenced), rep(1, 3), tolerance = 1e-4)
})

test_that("the stochastic runner, whose patch names no stiff block, refuses ARK", {
  f <- ark_fixture()
  expect_error(StochasticPatchRunner("TF24", "TF24_Env")(f$p, f$env, f$ctrl),
               "names a stiff block")
})
