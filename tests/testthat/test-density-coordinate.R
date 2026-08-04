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

## Multiple species share one competition profile: Patch::compute_competition()
## sums their per-species integrals into a single optical depth, each taken over
## that species' own node list. So the coordinate change has to hold for each
## species *through* a profile the others also shape. These use the established
## multi-species setups from test-strategy-ff16.R and test-strategy-k93.R.
test_that("multi-species runs converge across density coordinates", {
  cases <- list(
    FF16 = function() {
      p0 <- scm_base_parameters("FF16")
      add_strategies(p0, trait_matrix(c(0.0825, 0.2625), "lma"),
                     hyperpar = FF16_hyperpar,
                     birth_rate = list(11.99177, 16.51006))
    },
    K93 = function() {
      p0 <- scm_base_parameters("K93")
      p0$max_patch_lifetime <- 35.10667
      sp <- trait_matrix(c(0.042, 0.063, 0.052,
                           8.5e-3, 0.014, 0.015,
                           2.2e-4, 4.6e-4, 3e-4,
                           0.008, 0.008, 0.008,
                           1.8e-4, 4.4e-4, 5.1e-4,
                           1.4e-4, 2.5e-3, 8.8e-3,
                           0.044, 0.044, 0.044),
                         c("b_0", "b_1", "b_2", "c_0", "c_1", "d_0", "d_1"))
      add_strategies(p0, sp, birth_rate = c(20, 20, 20))
    }
  )

  for (x in names(cases)) {
    gaps <- lapply(0:2, function(r) {
      p <- interleave_schedule(cases[[x]](), r)
      h <- offspring_in_coordinate(p, x, FALSE)
      b <- offspring_in_coordinate(p, x, TRUE)
      abs(b - h) / abs(h)
    })

    ## Every species converges, not just the aggregate.
    expect_true(all(gaps[[2]] < gaps[[1]]))
    expect_true(all(gaps[[3]] < gaps[[2]]))
    expect_true(all(gaps[[1]] / gaps[[3]] > 4))
    expect_true(all(gaps[[3]] < 2e-3))
  }
})

## K93's three-species case largely excludes its first species -- its offspring
## production is ~90x below the dominant one. A marginal species is where a
## change to the shared profile could bite hardest, so check its *relative*
## error is no worse than the dominant species', rather than only that the
## totals agree.
test_that("a competitively excluded species is no worse off", {
  x <- "K93"
  p0 <- scm_base_parameters(x)
  p0$max_patch_lifetime <- 35.10667
  sp <- trait_matrix(c(0.042, 0.063, 0.052,
                       8.5e-3, 0.014, 0.015,
                       2.2e-4, 4.6e-4, 3e-4,
                       0.008, 0.008, 0.008,
                       1.8e-4, 4.4e-4, 5.1e-4,
                       1.4e-4, 2.5e-3, 8.8e-3,
                       0.044, 0.044, 0.044),
                     c("b_0", "b_1", "b_2", "c_0", "c_1", "d_0", "d_1"))
  p <- add_strategies(p0, sp, birth_rate = c(20, 20, 20))
  p <- interleave_schedule(p, 1)

  h <- offspring_in_coordinate(p, x, FALSE)
  b <- offspring_in_coordinate(p, x, TRUE)
  rel <- abs(b - h) / abs(h)

  ## Species 1 really is marginal, so the premise of the test holds.
  expect_lt(h[[1]] / max(h), 0.05)
  ## And its relative error is of the same order as the rest.
  expect_lt(rel[[1]], 3 * max(rel[-1]))
})

## Reporting boundary: whichever coordinate the solver carries internally, R is
## handed a density in *height*, so tidy_outputs.R's `density`,
## interpolate_to_heights() and the plots keep their meaning. The carried
## quantity is reported alongside rather than lost.
collected_run <- function(birth_date, lifetime = 20) {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- lifetime
  p <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- birth_date
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), ctrl)
  scm$collect <- TRUE
  scm$run()
  scm
}

## Compared in log space: |log N_b - log N_h| is the log of the density ratio,
## so 1e-3 means the densities agree to ~0.1%. A *relative* comparison of the logs
## is meaningless here because log_density crosses zero.
log_density_gap <- function(refine) {
  p_of <- function() {
    p0 <- scm_base_parameters("FF16")
    p0$max_patch_lifetime <- 20
    add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  }
  one <- function(birth_date) {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- birth_date
    scm <- SCM("FF16", "FF16_Env")(interleave_schedule(p_of(), refine),
                                   Environment("FF16"), ctrl)
    scm$run()
    scm$patch$species[[1]]
  }
  sa <- one(FALSE)
  sb <- one(TRUE)
  list(gap = abs(sb$log_densities - sa$log_densities),
       state_gap = abs(sb$log_densities_state - sa$log_densities))
}

test_that("log_densities is a height density in both coordinates", {
  g0 <- log_density_gap(0)
  g1 <- log_density_gap(1)

  ## The typical node agrees closely, and refining the schedule improves it at
  ## about the 2nd order of the differenced Jacobian (measured ratios ~3.9).
  expect_lt(median(g0$gap), 5e-3)
  expect_gt(median(g0$gap) / median(g1$gap), 2.5)
  expect_gt(quantile(g0$gap, 0.9) / quantile(g1$gap, 0.9), 2.5)

  ## Deliberately not asserted on max(): the *worst* node does not converge.
  ## |dh/dtau| is a ratio of two differences that both shrink as the schedule is
  ## refined, so cancellation error sets a floor, and refinement adds nodes in
  ## the near-empty tail (log density ~ -11) where that floor is worst. The
  ## reconstructed height density is sound in aggregate and for plotting; it
  ## should not be trusted node-by-node out in the tail.

  ## And the conversion is doing real work: unconverted, the carried quantity
  ## differs from the height density by ~1 in log space, not ~1e-3.
  expect_gt(median(g0$state_gap), 0.5)
})

test_that("log_densities_state is the quantity actually integrated", {
  a <- collected_run(FALSE)
  b <- collected_run(TRUE)
  sa <- a$patch$species[[1]]
  sb <- b$patch$species[[1]]

  ## On the height path there is nothing to convert.
  expect_equal(sa$log_densities_state, sa$log_densities)
  ## On the birth-date path it is the raw per-node state, which Node reports
  ## unconverted (a lone node cannot form dh/dtau).
  expect_equal(sb$log_densities_state,
               vapply(sb$nodes, function(nd) nd$log_density, numeric(1)))
  expect_false(isTRUE(all.equal(sb$log_densities_state, sb$log_densities)))

  ## And the two are related by the Jacobian, boundary node last.
  jac <- sb$height_jacobian
  expect_length(jac, sb$size + 1L)
  expect_equal(sb$log_densities,
               sb$log_densities_state - log(head(jac, sb$size)))
})

## The boundary node needs no differencing: dh/dtau = -g(H_0) at birth, and
## compute_initial_conditions() re-evaluates that every step, so the recorded
## rate is both exact and current for that node alone.
test_that("the boundary node's Jacobian is exact, not differenced", {
  b <- collected_run(TRUE)
  sb <- b$patch$species[[1]]
  n <- sb$size

  g0 <- sb$new_node$growth_rate_at_birth
  expect_gt(g0, 0)
  expect_identical(sb$height_jacobian[[n + 1L]], g0)

  ## And it is a real improvement over the one-sided difference it replaces:
  ## the boundary node sits a whole introduction interval from its neighbour.
  one_sided <- abs((sb$new_node$height - sb$heights[[n]]) /
                   (sb$new_node$introduction_time - sb$node_times[[n]]))
  expect_gt(abs(one_sided - g0) / g0, 0.1)

  ## Interior nodes are unaffected -- they never used the boundary node's own
  ## Jacobian, only its (h, tau) as the far point of a central difference.
  expect_true(all(is.finite(sb$height_jacobian[seq_len(n)])))
})

## The height coordinate's boundary condition is N(H_0) = birth_rate * pr_estab
## / g(H_0) (eq-bc1), and the birth-date one drops the division. Recording
## g(H_0) lets the two be checked against each other directly.
test_that("g(H_0) reproduces the height coordinate's boundary condition", {
  a <- collected_run(FALSE)
  sa <- a$patch$species[[1]]
  nd <- sa$new_node

  g0 <- nd$growth_rate_at_birth
  expect_gt(g0, 0)
  ## exp(log_density) * g == birth_rate * pr_estab, which the birth-date path
  ## carries directly. Compare in log space.
  expect_equal(nd$log_density + log(g0),
               log(1.0 * nd$individual$establishment_probability(a$patch$environment)),
               tolerance = 1e-10)
})

test_that("the collected state carries both densities", {
  a <- collected_run(FALSE)
  b <- collected_run(TRUE)

  rows_h <- rownames(a$patch$state$species[[1]])
  rows_b <- rownames(b$patch$state$species[[1]])
  expect_false("log_density_state" %in% rows_h)
  expect_true("log_density_state" %in% rows_b)

  ## The matrix that feeds tidy_patch() agrees with the accessor, so the
  ## `density = exp(log_density)` downstream is in height units.
  sb <- b$patch$species[[1]]
  st <- b$patch$state$species[[1]]
  expect_equal(unname(st["log_density", seq_len(sb$size)]), sb$log_densities)
  expect_equal(unname(st["log_density_state", seq_len(sb$size)]),
               sb$log_densities_state)
})

test_that("resume reads the raw state, so it is unaffected by the conversion", {
  x <- "FF16"; e <- "FF16_Env"
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- TRUE

  scm <- collected_run(TRUE)
  state <- export_patch_state(scm, step = max(2L, length(scm$history) %/% 2L))
  p2 <- set_initial_state(scm$parameters, state)

  scm2 <- SCM(x, e)(p2, Environment(x), ctrl)
  scm2$collect <- TRUE
  scm2$run()

  ## The seeded patch reproduces the exported ODE state exactly -- which it
  ## could not if the reporting conversion had leaked into the resume path.
  expect_equal(scm2$history[[1]]$ode_state, state$ode_state)
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
