# Carrying the SCM's size distribution in birth date rather than in height
# (Control$node_density_in_birth_date).
#
# The two coordinates describe the same population, so for a strategy whose
# growth is a function of size alone they must agree in the limit of a fine
# node schedule. They do not agree at a coarse one, and the gap is almost
# entirely the *height* coordinate's error: see the convergence test below,
# where the birth-date answer is already converged at the default schedule
# while the height answer is still climbing toward it.

interleave_schedule <- function(p, n = 1) {
  for (i in seq_len(n)) {
    p$node_schedule_times <- lapply(p$node_schedule_times, function(tt) {
      sort(unique(c(tt, (head(tt, -1) + tail(tt, -1)) / 2)))
    })
  }
  p
}

size_only_parameters <- function(x) {
  p0 <- scm_base_parameters(x)
  if (x == "FF16") {
    add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  } else {
    add_strategies(p0, trait_matrix(0.059, "b_0"), birth_rate = 1.0)
  }
}

offspring_in_coordinate <- function(p, x, birth_date, step_max = NULL) {
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- birth_date
  if (!is.null(step_max)) {
    ctrl$ode_step_size_max <- step_max
  }
  run_scm(p, Environment(x), ctrl)$offspring_production
}

## The load-bearing claim: for FF16 and K93, growth is a function of size, so
## the compression term the height coordinate computes is correct and the two
## coordinates must converge to the same offspring production. A stale or
## duplicated birth-date abscissa breaks this, which is what makes it a useful
## regression test rather than a restatement of the implementation.
test_that("size-only strategies converge across density coordinates", {
  for (x in c("FF16", "K93")) {
    p <- size_only_parameters(x)

    gaps <- vapply(0:2, function(r) {
      pr <- interleave_schedule(p, r)
      h <- offspring_in_coordinate(pr, x, FALSE)
      b <- offspring_in_coordinate(pr, x, TRUE)
      abs(b - h) / abs(h)
    }, numeric(1))

    ## Halving the node spacing must shrink the gap, and at roughly the
    ## second order of the trapezium rule (measured ratios 2.4, 3.6, 4.0 on
    ## FF16). Loose bounds: the point is convergence, not the exact rate.
    expect_true(all(diff(gaps) < 0))
    expect_gt(gaps[[1]] / gaps[[3]], 4)

    ## And they really are converging to each other, not merely narrowing.
    expect_lt(gaps[[3]], 2e-3)
  }
})

## The birth-date answer is the accurate one: it is what both coordinates
## converge to, and it gets there on a far coarser schedule. This is the
## substantive claim of the change, so pin it rather than leaving it in prose.
test_that("the birth-date coordinate converges on a coarser schedule", {
  x <- "FF16"
  p <- size_only_parameters(x)

  refined <- offspring_in_coordinate(interleave_schedule(p, 3), x, TRUE)
  coarse_birth <- offspring_in_coordinate(p, x, TRUE)
  coarse_height <- offspring_in_coordinate(p, x, FALSE)

  ## At the default schedule the height coordinate is ~50x further from the
  ## converged value than the birth-date coordinate is.
  err_birth <- abs(coarse_birth - refined) / refined
  err_height <- abs(coarse_height - refined) / refined
  expect_lt(err_birth, 1e-3)
  expect_gt(err_height / err_birth, 10)
})

## The boundary node's birth date is the current time, and compute_rates()
## (which stamps it) runs after the set_ode_state() that rebuilds the
## environment. Reading the stamp during that rebuild would use the previous
## derivs call's time, making the last trapezium interval a function of the
## step size. The measured effect on offspring production is small (<1e-6 on
## FF16), so assert the invariant directly rather than through a result.
test_that("the boundary node's birth date tracks patch time", {
  for (x in c("FF16", "K93")) {
    e <- switch(x, FF16 = "FF16_Env", K93 = "K93_Env")
    p <- size_only_parameters(x)
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- TRUE
    patch <- Patch(x, e)(p, Environment(x), ctrl)

    for (t in c(0.0, 1.5, 7.25)) {
      patch$set_time(t)
      patch$compute_environment()
      expect_identical(patch$species[[1]]$new_node$introduction_time, t)
    }
  }
})

## Introduction times are the birth-date coordinate's quadrature grid, so
## repeated ones span zero width and drop silently out of the integrals. A
## scheduled run cannot produce them; a patch seeded without per-node times
## can, and used to do so silently.
test_that("nodes sharing a birth date are rejected in birth-date coordinates", {
  x <- "FF16"; e <- "FF16_Env"
  p <- size_only_parameters(x)

  scm <- SCM(x, e)(p, Environment(x), Control())
  scm$collect <- TRUE
  scm$run()
  state <- export_patch_state(scm, step = max(2L, length(scm$history) %/% 2L))

  ## A faithful resume carries per-node times, so both coordinates accept it.
  p_ok <- set_initial_state(p, state)
  ctrl_bd <- Control()
  ctrl_bd$node_density_in_birth_date <- TRUE
  expect_no_error(Patch(x, e)(p_ok, Environment(x), ctrl_bd))

  ## Drop them and every node inherits the boundary node's birth date.
  p_bad <- p_ok
  p_bad$initial_node_times <- numeric(0)
  expect_error(Patch(x, e)(p_bad, Environment(x), ctrl_bd),
               "sharing an introduction time")

  ## The height coordinate never integrates over these, so it is unaffected.
  expect_no_error(Patch(x, e)(p_bad, Environment(x), Control()))
})

## The coordinate is stored per strategy but cannot differ between the species
## of one patch: Patch::add_strategies() overwrites every strategy's control
## with the patch's. Pin that, because compute_competition() sums the species'
## contributions into a single optical depth and would otherwise be adding a
## birth-date integral to a height integral.
test_that("the patch control decides the coordinate for every species", {
  x <- "FF16"; e <- "FF16_Env"
  ctrl_bd <- Control()
  ctrl_bd$node_density_in_birth_date <- TRUE

  p0 <- scm_base_parameters(x)
  p <- add_strategies(p0, trait_matrix(c(0.08, 0.2), "lma"),
                      birth_rate = list(1.0, 1.0))
  ## Deliberately disagree with the patch control; it should be ignored.
  p$strategies[[1]]$control <- Control()

  patch <- Patch(x, e)(p, Environment(x), ctrl_bd)
  for (s in patch$species) {
    expect_true(s$density_in_birth_date)
  }

  ## And the other way round: a patch built with the default control carries
  ## every species in height, whatever the strategies say.
  p$strategies[[1]]$control <- ctrl_bd
  patch_h <- Patch(x, e)(p, Environment(x), Control())
  for (s in patch_h$species) {
    expect_false(s$density_in_birth_date)
  }
})
