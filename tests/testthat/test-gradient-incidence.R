# The instruments the gradient's scope questions are answered with.
#
# The leaf classifies its operating point by the branch taken and the next plant
# overwrites it, and a clamp that severs a row leaves a number indistinguishable
# from a true zero. Neither is recoverable after a run from anything else, so
# these counters are the only route to "how often" -- and "how often" is what
# decides whether a refused regime is a corner or most of the run.

incidence_stand <- function(rain, lifetime, k_I = 0.5) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  tr <- c(lma = 0.0825, hmat = 5.13, k_I = k_I, a_l1 = 5.44, a_l2 = 0.306)
  p <- add_strategies(p, trait_matrix(unname(tr), names(tr)),
                      hyperpar = TF24_hyperpar, birth_rate = list(1.10))
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, 5))
  env$extrinsic_drivers_set_constant("rainfall", rain)
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- TRUE
  scm <- SCM("TF24", "TF24_Env")(p, env, ctrl)
  census_clear_operating_point_counts_tf24(scm)
  scm$run()
  scm
}

incidence_of <- function(scm) {
  counts <- census_operating_point_counts_tf24(scm)[[1]]
  stats::setNames(counts, census_operating_point_names_tf24())
}

test_that("the classification tally is the route to a regime's incidence", {
  # A wet stand never leaves the branch the gradient answers for, which is what
  # makes it the fixture every other rung uses -- and is why incidence measured
  # on one says nothing about a dry one.
  wet <- incidence_of(incidence_stand(2.0, 5))
  expect_gt(wet[["interior"]], 0)
  expect_equal(sum(wet[names(wet) != "interior"]), 0)

  # The tally is cleared and re-accumulated per run rather than carried, or a
  # second measurement would read the first one's states as well.
  scm <- incidence_stand(2.0, 5)
  census_clear_operating_point_counts_tf24(scm)
  expect_equal(sum(incidence_of(scm)), 0)
})

test_that("a refused gradient is a small minority of the operating points", {
  # The number this exists to produce. A refusal is metric-level and total, so
  # the share of points that cause it is what says how much answering the pinned
  # branch would buy -- and it is not recoverable from the refusal message, which
  # names the first such point and nothing about the rest.
  scm <- incidence_stand(0.25, 10)
  n <- incidence_of(scm)
  total <- sum(n)
  dry <- n[["pinned-dry-root-crit"]] + n[["pinned-dry-root-psi-crit"]]
  expect_gt(dry, 0)
  expect_gt(n[["interior"]], dry)

  # Which arm bound is not a detail: the two are different functions of the
  # inputs, so the row a pinned point needs depends on it. At shipped defaults
  # the root's own critical potential never wins the min, which is why a fixture
  # for that arm has to lower it deliberately rather than wait for one.
  expect_equal(n[["pinned-dry-root-psi-crit"]], 0)
  expect_gt(n[["pinned-dry-root-crit"]], 0)

  share <- 100 * dry / total
  message(sprintf("  dry pins on a refusing run: %.0f of %.0f solves (%.2f%%), all on the %s arm",
                  dry, total, share, "root-crit"))
  # Stated as a band rather than a value: the exact count moves with the
  # schedule the adaptive pass resolves, and what is being pinned is that a
  # total refusal is caused by a small minority of points.
  expect_lt(share, 5)

  # And the pairing that makes the number mean something: this run's gradient
  # really is refused, so the minority above is what costs the whole answer.
  g <- stand_gradient(scm)
  expect_true(all(g$status == "refused"))
  expect_false(is.null(g$refusal[[1]]))
})

test_that("the light floor is counted, and does not bind at shipped values", {
  # Where this clamp binds, a cohort's radiation stops depending on any other
  # cohort's height and the row is severed by the guard rather than by the model.
  # Both halves are asserted, because a counter that never fires and a counter
  # that always fires are equally uninformative.
  shipped <- census_clamp_counts_tf24(incidence_stand(2.0, 5, k_I = 0.5))[[1]]
  expect_equal(shipped[[1]], 0)

  # k_I is a free parameter a gradient-driven search walks, and walking it up is
  # what walks the field into the floor. Eighty times the shipped value.
  walked <- incidence_stand(2.0, 5, k_I = 40)
  fired <- census_clamp_counts_tf24(walked)[[1]][[1]]
  solves <- sum(incidence_of(walked))
  message(sprintf("  light floor at k_I = 40: %.0f of %.0f solves (%.1f%%)",
                  fired, solves, 100 * fired / solves))
  expect_gt(fired, 0)

  expect_identical(census_clamp_names_tf24(), "light_floor")
})
