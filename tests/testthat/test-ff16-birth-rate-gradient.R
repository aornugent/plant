# AD-10 (#15): the FF16 birth-rate axis on the coupled resident replay. A single
# registered initial-state leaf -- the focal species' birth-rate scale -- is seeded
# active; it multiplies the extrinsic birth rate at every node introduction, so it
# flows through the cohort density and, because the canopy is recomputed live on the
# recorded knots, through the whole stand. The reverse-mode derivative must match a
# central finite difference of the same frozen-knot replay (perturb the scale, read
# the metric value). Exercises the compiled birth_rate_gradient path in plant.so.
#
# The birth-rate axis differs from the trait/invasion axes in one numerical respect:
# the coupled replay rebuilds the canopy through the SCM's data-dependent trapezium
# competition integral (its early-exit branches), so the replay's output as a function
# of the seed is only piecewise smooth -- isolated O(1e-5) kinks. The reverse-mode AD
# is the exact derivative of that piecewise-smooth map; central FD converges to it as
# h shrinks (below), but the census integrals over the whole active size distribution
# amplify the kinks to a ~1e-3 FD noise floor, so they are gated more loosely than the
# scalar offspring / R0 metrics. The invasion (frozen-canopy) axis has no such kink.

test_that("FF16 birth-rate gradient matches central FD of the coupled resident replay", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE

  b <- 2  # non-unit base rate exercises the d/d(scale) -> d/d(birth_rate) rescale
  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = b)
  scm <- run_scm(pr1, e, ctrl)

  metrics <- c("offspring_production", "net_reproduction_ratio",
               "biomass", "LAI", "basal_area")

  ad <- birth_rate_gradient_cpp(scm, metrics, 1L, "FF16", 1.0)
  expect_equal(ad$birth_rate, b)
  expect_true(all(is.finite(ad$gradient)))
  expect_named(ad$gradient, metrics)

  # Central FD of the SAME frozen-knot replay: the metric value at birth-rate scale
  # s is the metric at effective rate s * b, so a central difference over the scale
  # (divided by b) is d(metric)/d(birth_rate). Same replay both sides, so this
  # isolates the AD's exactness from the record->replay residual.
  h <- 1e-5
  vplus <- birth_rate_gradient_cpp(scm, metrics, 1L, "FF16", 1 + h)$value
  vminus <- birth_rate_gradient_cpp(scm, metrics, 1L, "FF16", 1 - h)$value
  fd <- (vplus - vminus) / (2 * h * b)

  # The scalar birth-rate deliverables (the offspring axis and dR0/db) match to the
  # AD-vs-FD tolerance; the census metrics carry the SCM-kink noise floor.
  for (mt in c("offspring_production", "net_reproduction_ratio")) {
    expect_equal(unname(ad$gradient[mt]), unname(fd[mt]), tolerance = 1.5e-4)
  }
  for (mt in c("biomass", "LAI", "basal_area")) {
    expect_equal(unname(ad$gradient[mt]), unname(fd[mt]), tolerance = 1e-3)
  }
})

test_that("the resident birth-rate feedback flips the sign against the frozen identity", {
  # The frozen part of the birth-rate response is the identity metric / birth_rate
  # (density scales linearly with the seed rain). The coupled canopy feedback is
  # non-trivial: for offspring production it overwhelms and flips the sign -- more
  # seed rain shades the stand enough that total survival-weighted offspring falls.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control(); ctrl$save_RK45_cache <- TRUE
  b <- 2
  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = b)
  scm <- run_scm(pr1, e, ctrl)

  g <- birth_rate_gradient(scm, "offspring_production")
  identity <- g$value[["offspring_production"]] / b   # frozen-canopy d/db
  expect_gt(identity, 0)                                # naive scaling is positive
  expect_lt(g$gradient[["offspring_production"]], 0)    # feedback flips it negative
})

test_that("a Newton step using dR0/db moves the resident toward R0 = 1", {
  # dR0/db is the plant-side derivative for the equilibrium (R0 = 1) Newton solve.
  # One step b - (R0 - 1) / (dR0/db) from an off-equilibrium birth rate must reduce
  # |R0 - 1| on a re-run resident (correct sign and magnitude of the density feedback).
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control(); ctrl$save_RK45_cache <- TRUE
  b <- 2

  run_R0 <- function(birth) {
    pr <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = birth)
    run_scm(pr, e, ctrl)$net_reproduction_ratios[1]
  }

  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = b)
  scm <- run_scm(pr1, e, ctrl)
  R0 <- scm$net_reproduction_ratios[1]
  dR0db <- birth_rate_gradient(scm, "net_reproduction_ratio")$gradient[["net_reproduction_ratio"]]

  expect_lt(dR0db, 0)                 # more seed rain -> more shading -> lower R0
  expect_false(isTRUE(all.equal(R0, 1, tolerance = 1e-2)))  # genuinely off-equilibrium

  b_new <- b - (R0 - 1) / dR0db
  expect_gt(b_new, 0)
  R0_new <- run_R0(b_new)
  expect_lt(abs(R0_new - 1), abs(R0 - 1))
})

test_that("multi-species birth-rate gradient is finite and matches FD for the focal species", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control(); ctrl$save_RK45_cache <- TRUE

  b <- c(2, 1.5)
  pr <- add_strategies(p0, trait_matrix(c(0.1, 0.2), "lma"), birth_rate = b)
  scm <- run_scm(pr, e, ctrl)

  mets <- c("offspring_production", "biomass")
  ad <- birth_rate_gradient_cpp(scm, mets, 2L, "FF16", 1.0)  # differentiate species 2
  expect_equal(ad$birth_rate, b[2])
  expect_true(all(is.finite(ad$gradient)))

  # Seeding only species 2's birth rate re-shades the whole (two-species) stand
  # through its density; the biomass derivative carries that coupled feedback and
  # matches FD cleanly. Offspring production for species 2 is a near-cancellation of
  # the identity by the feedback (a small residual), so it is checked in absolute
  # terms rather than relative -- the AD and FD agree to the census-kink floor.
  h <- 3e-5
  vplus <- birth_rate_gradient_cpp(scm, mets, 2L, "FF16", 1 + h)$value
  vminus <- birth_rate_gradient_cpp(scm, mets, 2L, "FF16", 1 - h)$value
  fd <- (vplus - vminus) / (2 * h * b[2])
  expect_equal(unname(ad$gradient[["biomass"]]), unname(fd[["biomass"]]), tolerance = 1e-3)
  expect_lt(abs(ad$gradient[["offspring_production"]] - fd[["offspring_production"]]), 1e-3)
})
