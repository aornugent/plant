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

test_that("the dry pins are a small minority, and the run answers over them", {
  # The number this exists to produce. It was first taken while this driver's
  # gradient was refused outright, to say how much answering the pinned branch
  # would buy; the branch answers now, so the same number says what the answer
  # rests on. It is not recoverable from a run afterwards either way -- the
  # classification is overwritten by the next individual.
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
  message(sprintf("  dry pins: %.0f of %.0f solves (%.2f%%), all on the %s arm",
                  dry, total, share, "root-crit"))
  # Stated as a band rather than a value: the exact count moves with the schedule
  # the adaptive pass resolves.
  expect_lt(share, 5)

  # And the pairing that makes the number mean something. This assertion used to
  # run the other way -- the run refused, and the minority above was what cost the
  # whole answer. A driver that reaches the pinned branch is now the driver that
  # demonstrates it, so what is asserted is that it comes back answered with
  # finite rows rather than that it comes back refused.
  g <- stand_gradient(scm)
  expect_true(all(g$status != "refused"))
  expect_null(g$refusal[[1]])
  expect_true(any(is.finite(g$gradient[[1]])))
})

test_that("the light floor is counted on both paths, and binds at neither shipped value", {
  # Where this clamp binds, a cohort's radiation stops depending on any other
  # cohort's height. Both halves are asserted, because a counter that never fires
  # and a counter that always fires are equally uninformative.
  #
  # ⚠️ THE FLOOR IS TWO SITES, NOT ONE, and the uncounted one bound first. The
  # crown site floors each quadrature point of the mean-light integrand; the
  # radiation site floors the mean those points make. Since the shape integrates
  # to one, a floored point cannot pull the mean below the floor -- so on the
  # shipped shading model the crown site fires ~150x more often, and reading the
  # radiation site alone under-reports the severance by that factor.
  # Scoped to the two light sites: the other clamps have their own regimes and
  # their own block below, and several of them do bind at shipped values.
  nm <- census_clamp_names_tf24()
  light <- match(c("light_floor", "light_floor_crown"), nm)
  shipped <- census_clamp_counts_tf24(incidence_stand(2.0, 5, k_I = 0.5))[[1]]
  expect_true(all(shipped[light] == 0))

  # k_I is a free parameter a gradient-driven search walks, and walking it up is
  # what walks the field into the floor. Eighty times the shipped value.
  walked <- incidence_stand(2.0, 5, k_I = 40)
  fired <- census_clamp_counts_tf24(walked)[[1]]
  solves <- sum(incidence_of(walked))
  message(sprintf("  at k_I = 40, over %.0f solves: %s", solves,
                  paste(sprintf("%s %.0f (%.1f%%)", nm[light], fired[light],
                                100 * fired[light] / solves), collapse = "  ")))
  expect_true(all(fired[light] > 0))
  # The crown site binds first, so it cannot be the smaller of the two.
  expect_gt(fired[[light[[2]]]], fired[[light[[1]]]])

  # The forward model keeps running AND the gradient answers: below the floor the
  # census is not a function of light at all, so the row is exactly zero for the
  # model as evaluated rather than a row withheld.
  expect_gt(stand_census(walked)[[1]], 0)
  g <- stand_gradient(walked)
  expect_true(all(g$status != "refused"))
  expect_true(all(is.finite(g$gradient)))

  # And the severance is readable rather than silent, which is the whole basis on
  # which the zero is declared instead of refused. The forward tally cannot stand
  # in for this: it counts every solve, where the sweep visits only the recorded
  # steps.
  swept <- census_clamp_counts_differentiated_tf24(walked)[[1]]
  message(sprintf("  the sweep's own severances: %s",
                  paste(sprintf("%s %.0f", nm[light], swept[light]),
                        collapse = "  ")))
  expect_true(all(swept[light] > 0))
  # A counter a rebind drops reports zero however often the clamp fires, which is
  # what this one did before the storage was shared. And the sweep visits the
  # recorded steps rather than every solve, so its tally is the smaller.
  expect_lt(swept[[light[[1]]]], fired[[light[[1]]]])
  expect_lt(swept[[light[[2]]]], fired[[light[[2]]]])

  # Named from the enum, so a site cannot be counted under its neighbour's name.
  expect_true(all(c("light_floor", "light_floor_crown") %in% nm))
  expect_length(fired, length(nm))
  expect_length(swept, length(nm))
})

# Every clamp site, classified by what its incidence says rather than by what
# reading the code suggests. The classification lives here rather than in a
# document because a document cannot fail.
#
# ⚠️ THE TEST IS THREE-WAY AND THE FIRST CASE SPLITS IN TWO, which reading the
# code does not reveal:
#
#   the model's own zero   the clamp IS a modelling statement, so the census
#                          genuinely does not depend on what it masks. Declared,
#                          and NOT a candidate for removal -- rooting_depth is
#                          this, and it is the largest severance in the model.
#   a guard's zero         a numerical floor, and the census is bit-identical
#                          either side of it, so the row is exactly zero for the
#                          model as evaluated. Declared, and a candidate for
#                          removal by changing the FORWARD model.
#   never binds            counted, and reported as never having fired, which is
#                          the only thing that separates a guard that held from
#                          one nothing reached.
clamp_class <- list(
  rooting_depth          = "model",
  light_floor            = "guard",
  light_floor_crown      = "guard",
  storage_floor          = "guard",
  reserve_ceiling        = "guard",
  soil_moisture_floor    = "forward-only",
  soil_potential_ceiling = "forward-only",
  soil_conductivity      = "forward-only",
  soil_positivity        = "never",
  rainfall               = "never",
  infiltration           = "never"
)

test_that("every clamp site is classified, and by a measured incidence", {
  nm <- census_clamp_names_tf24()
  # A site with no classification is the drift this list exists to prevent.
  expect_setequal(nm, names(clamp_class))

  swept_of <- function(scm) {
    stats::setNames(census_clamp_counts_differentiated_tf24(scm)[[1]], nm)
  }
  fwd_of <- function(scm) {
    stats::setNames(census_clamp_counts_tf24(scm)[[1]], nm)
  }

  wet <- incidence_stand(2.0, 5)
  invisible(stand_gradient(wet))
  wet_s <- swept_of(wet)
  dry <- incidence_stand(0.10, 5)
  invisible(stand_gradient(dry))
  dry_s <- swept_of(dry)
  dry_f <- fwd_of(dry)

  message(sprintf("  %-24s %10s %10s", "site", "wet", "drought"))
  for (s in nm) {
    message(sprintf("  %-24s %10.0f %10.0f  [%s]", s, wet_s[[s]], dry_s[[s]],
                    clamp_class[[s]]))
  }

  # Never means never, on both paths and both drivers. A site in this class that
  # starts firing is a regime the suite has never seen.
  for (s in names(clamp_class)[unlist(clamp_class) == "never"]) {
    expect_equal(wet_s[[s]], 0, label = paste(s, "fired on the wet driver"))
    expect_equal(dry_s[[s]], 0, label = paste(s, "fired on the drought driver"))
    expect_equal(dry_f[[s]], 0, label = paste(s, "fired forward on drought"))
  }

  # Forward-only means the sweep never met it, which is a weaker claim than
  # never and has to be kept separate: the guard is reachable, so a longer run
  # or a finer schedule could put it on a recorded step.
  for (s in names(clamp_class)[unlist(clamp_class) == "forward-only"]) {
    expect_equal(dry_s[[s]], 0,
                 label = paste(s, "reached the differentiated path"))
  }

  # The model's own zero binds on EVERY driver including the control, which is
  # what says it is the model rather than a corner.
  expect_gt(wet_s[["rooting_depth"]], 0)
  expect_gt(dry_s[["rooting_depth"]], 0)

  # And a guard's zero binds where its own regime is reached and not on the
  # control, which is what makes it a guard rather than the model.
  expect_gt(dry_s[["storage_floor"]], 0)
  expect_equal(wet_s[["storage_floor"]], 0)
  expect_equal(wet_s[["light_floor"]], 0)
  expect_equal(wet_s[["light_floor_crown"]], 0)

  # Both drivers still answer, carrying those declared zeros. A severance that
  # made the gradient wrong would have to show up as a refusal or as a
  # disagreement with a rebuilt difference, and neither is here.
  expect_true(all(as.vector(stand_gradient(wet)$status) != "refused"))
  expect_true(all(as.vector(stand_gradient(dry)$status) != "refused"))
})

test_that("phylloptim's root-vulnerability clamps stay out of reach, with a margin", {
  # These are NOT counted, deliberately, and this is what stands in for a counter.
  #
  # Two clamps in the root curve are thresholds on a layer potential: the curve's
  # argument is clamped into its knot domain at 6.8229 MPa, and the cumulative
  # integral is capped past 7.3132. Both already carry matching derivative kills,
  # so where they bind the row is an honest zero rather than a wrong number -- and
  # the root_b row stays right under the cap by the homogeneity rather than in
  # spite of it. So what is worth asserting is not a count but the DISTANCE, which
  # a counter reading zero cannot report.
  #
  # Instrumenting them would need a phylloptim header edit, hence a reinstall,
  # hence a near-full plant recompile: the cost is the build loop rather than the
  # code. A threshold on a readable quantity does not need it.
  last_knot <- 6.8229
  integral_cap <- 7.3132

  worst <- 0
  for (d in list(list(rain = 2.00, name = "wet"),
                 list(rain = 0.10, name = "drought"),
                 list(rain = 0.05, name = "very-dry"))) {
    scm <- incidence_stand(d$rain, 5)
    e <- scm$patch$environment
    nlayer <- e$get_soil_number_of_depths()
    # Over the RECORDED steps, which is the set the sweep visits -- a terminal
    # reading misses a layer that dried and rewetted, and those are the states the
    # clamp would bind in.
    for (r in scm$store_trajectory()) {
      s <- r$state
      n <- length(s)
      theta <- s[(n - 3 - nlayer):(n - 4)]
      if (length(theta) != nlayer || any(!is.finite(theta))) next
      psi <- vapply(theta, function(x) e$psi_from_soil_moist(x), numeric(1))
      if (any(!is.finite(psi))) next
      worst <- max(worst, max(psi))
    }
  }
  message(sprintf("  worst layer potential over three drivers: %.4f MPa, against %.4f and %.4f",
                  worst, last_knot, integral_cap))
  expect_lt(worst, last_knot)
  expect_lt(worst, integral_cap)

  # Non-vacuity: a run that dried nothing would pass the two bounds above while
  # measuring nothing at all, so the drying has to be real.
  expect_gt(worst, 1.0)
})

test_that("the curvature guard reports how close it came, not only that it held", {
  # A guard that held and a guard nothing reached report the same green, so the
  # distance to the floor is carried out of the run. The floor is a declared
  # Control entry for the same reason: it changes which rows exist, so two
  # gradients taken at different values are gradients of different functions.
  wet <- incidence_stand(2.0, 5)
  invisible(stand_gradient(wet))
  margin <- census_curvature_margin_tf24(wet)[[1]]
  floor <- gradient_control(wet)[["gradient_curvature_floor"]]
  message(sprintf("  smallest curvature met: %.4g, against a floor of %.4g (%.0fx)",
                  margin, floor, margin / floor))
  # Non-vacuity both ways: a margin of -1 means the interior branch was never
  # reached, so the reading would say nothing.
  expect_gt(margin, 0)
  expect_gt(margin, floor)
  # And the floor is in the set stand_gradient compares, or two gradients taken
  # at different floors would read as comparable.
  expect_true("gradient_curvature_floor" %in% names(gradient_control(wet)))
})
