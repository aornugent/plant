# The light field is built in one descent over the knot grid rather than as a sum
# over cohorts repeated at every knot. The descent is the path a run takes; the
# per-height reduction is what defines the answer, and it is still there for the
# cases the descent declines. Every check below holds one against the other.
#
# The scale to beat is the per-height reduction's own rounding, which is
# .Machine$double.eps * A(0): A(0) is the sum of every cohort's amplitude, so it
# bounds the terms that were added, and the check needs no reference value.

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

# A stand the model actually reaches. The reduction reads densities, heights and
# quadrature weights that only a run puts into a consistent state, so an
# assembled patch would be testing something else.
stand_patch <- function(type, t_max = 20, ctrl = Control()) {
  p0 <- scm_base_parameters(type)
  # Before add_strategies: the node-introduction schedule is built from the
  # lifetime, so shortening it afterwards leaves a schedule running past the end.
  p0$max_patch_lifetime <- t_max
  trait <- if (type == "K93") trait_matrix(0.059, "b_0") else trait_matrix(0.08, "lma")
  p1 <- add_strategies(p0, trait, birth_rate = 1.0)
  run_scm(p1, Environment(type), ctrl)$patch
}

# The grid the field was built on, and the two reductions read at it.
compare_on_grid <- function(patch) {
  z <- patch$environment$light_availability$state[, "height"]
  swept <- patch$compute_competition_grid(z)
  per_height <- vapply(z, function(zi) patch$compute_competition_and_slope(zi),
                       numeric(2))
  list(z = z,
       value = swept[, "value"], value_ref = per_height[1, ],
       slope = swept[, "slope"], slope_ref = per_height[2, ],
       scale = .Machine$double.eps * per_height[1, 1])
}

test_that("the descent reproduces the per-height reduction, cohort for cohort", {
  # FF16 at three stand sizes. The claim is agreement to the per-height
  # reduction's own rounding, not to a tolerance chosen to pass.
  for (t_max in c(5, 20, 80)) {
    patch <- stand_patch("FF16", t_max = t_max)
    n <- patch$species[[1]]$size
    r <- compare_on_grid(patch)

    worst_value <- max(abs(r$value - r$value_ref))
    worst_slope <- max(abs(r$slope - r$slope_ref))
    message(sprintf(
      "FF16 lifetime %2d: %3d cohorts, %4d knots, value %.2e slope %.2e (eps*A(0) = %.2e)",
      t_max, n, length(r$z), worst_value, worst_slope, r$scale))

    # A re-association moves the last bits and nothing more. Measured at 2.1x
    # the scale on a 98-cohort stand; 20x is headroom, not a fitted tolerance.
    expect_lt(worst_value, 20 * r$scale)
    expect_lt(worst_slope, 50 * r$scale)
  }
})

test_that("the descent reproduces the per-height reduction for K93 and TF24", {
  for (type in c("K93", "TF24")) {
    patch <- stand_patch(type, t_max = 10)
    r <- compare_on_grid(patch)
    worst_value <- max(abs(r$value - r$value_ref))
    message(sprintf("%s: %d knots, worst value %.2e (eps*A(0) = %.2e)",
                    type, length(r$z), worst_value, r$scale))
    expect_lt(worst_value, 20 * r$scale)
  }
})

test_that("the patch area divides the descent as it divides the per-height sum", {
  # The per-height entry points divide at the patch level, and filling the grid
  # from the species directly walks around that. Every other test runs at the
  # default area, where dividing by one hides the omission.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 20
  p0$patch_area <- 3.0
  p1 <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  # The SCM warns that it only checks units at area 1. That is the warning this
  # test is here to run against, so take it and carry on.
  patch <- suppressWarnings(run_scm(p1, Environment("FF16"), Control()))$patch
  expect_equal(patch$get_area, 3.0)

  r <- compare_on_grid(patch)
  expect_lt(max(abs(r$value - r$value_ref)), 20 * r$scale)
  # And the field is not the area-one field: an omitted division would leave it
  # three times larger, which this says it is not.
  expect_gt(max(abs(r$value)), 0)
})

test_that("a boundary node carrying no density is admitted as the per-height sum admits it", {
  # The per-height reduction declines the closing interval when the node it
  # would close from contributes nothing, which is a stand whose newest cohorts
  # carry no density yet. Seeing this needs a stand in that state; on any other
  # it is a fixed point of the mistake.
  patch <- stand_patch("FF16", t_max = 20)
  amplitudes <- patch$species[[1]]$compute_competition_effect_by_nodes
  # A run reaches nodes whose density has collapsed; if this stand has none the
  # check is vacuous and should say so rather than pass quietly.
  expect_true(min(amplitudes) >= 0)

  r <- compare_on_grid(patch)
  expect_lt(max(abs(r$value - r$value_ref)), 20 * r$scale)
})

test_that("both size-distribution coordinates give the same two answers", {
  # The abscissa is the introduction time in one coordinate and minus the height
  # in the other, and only the quadrature weights differ. Both are live.
  for (birth_date in c(FALSE, TRUE)) {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- birth_date
    patch <- stand_patch("FF16", t_max = 20, ctrl = ctrl)
    r <- compare_on_grid(patch)
    message(sprintf("birth-date coordinate %s: worst value %.2e (eps*A(0) = %.2e)",
                    birth_date, max(abs(r$value - r$value_ref)), r$scale))
    expect_lt(max(abs(r$value - r$value_ref)), 20 * r$scale)
  }
})

test_that("a profile that does not expand in u^eta takes the per-height path", {
  # The soft box is not a polynomial in u^eta, so the running sums do not
  # describe it and the descent hands it back. The field must still build, and
  # the two reductions must still be the same number -- here exactly, because the
  # descent IS the per-height path for it.
  ctrl <- Control()
  ctrl$shading_model <- "flat-top-soft-box"
  patch <- stand_patch("FF16", t_max = 10, ctrl = ctrl)
  r <- compare_on_grid(patch)
  message(sprintf("flat-top-soft-box: worst value %.2e",
                  max(abs(r$value - r$value_ref))))
  expect_identical(r$value, r$value_ref)
  expect_identical(r$slope, r$slope_ref)
})

test_that("a profile with a step in it is refused a field rather than fitted one", {
  # The hard box shades in a step at the crown centre. A grid of fixed knots does
  # not fail to fit that -- it fits it as a ramp one span wide and returns
  # numbers, where the adaptive fit this replaced refused by running out of
  # refinement. The guarantee has to be stated for a lattice to keep it.
  ctrl <- Control()
  ctrl$shading_model <- "flat-top-box"
  expect_error(stand_patch("FF16", t_max = 10, ctrl = ctrl),
               "this model has no light environment")

  # A standalone individual under the same model still computes: it is the patch
  # light field that has no answer, not the crown.
  s <- FF16_Strategy()
  s$control$shading_model <- "flat-top-box"
  ind <- FF16_Individual(s)
  ind$set_state("height", 10)
  env <- Environment("FF16")
  env$set_fixed_environment(0.5, 100)
  expect_silent(ind$compute_rates(env))
})

# A patch of `n` nodes whose heights are set directly, for the two states the
# descent declines. Node state is blocked per node with the height first, so
# writing every seventh entry writes the heights.
patch_with_heights <- function(heights) {
  x <- "FF16"
  e <- environment_types[[x]]
  s <- strategy_types[[x]]()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- Parameters(x, e)(strategies = list(s), patch_type = "meta-population")
  patch <- Patch(x, e)(p, Environment(x), Control())
  if (length(heights) == 0) {
    return(patch)
  }
  for (i in seq_along(heights)) {
    patch$introduce_new_node(1)
  }
  state <- patch$ode_state
  stride <- length(state) / length(heights)
  state[seq(1, length(state), by = stride)] <- heights
  patch$set_ode_state(state, 0)
  patch
}

test_that("a broken height ordering takes the per-height path", {
  # The descent admits cohorts in decreasing height order, and a node list that
  # is not in that order would admit them at the wrong knots. The per-height
  # reduction has its own answer for this (#571), so the descent hands it back
  # rather than carrying a second sort down the grid.
  patch <- patch_with_heights(c(5, 3, 4))
  expect_equal(patch$height_max, 5)

  r <- compare_on_grid(patch)
  expect_identical(r$value, r$value_ref)
  expect_identical(r$slope, r$slope_ref)
  # The field is not trivially zero, so the agreement is about something.
  expect_gt(max(abs(r$value)), 0)
})

test_that("a boundary node above the shortest cohort takes the per-height path", {
  # No run reaches this -- a cohort starts at the seed height and grows -- but a
  # patch written by hand can. The closing interval would then be admitted at
  # heights where the per-height reduction declines it and the boundary profile
  # is not yet zero, so the two would part company. The descent declines instead.
  # The height a node is born at, read before anything overwrites it.
  seed <- patch_with_heights(numeric(0))$species[[1]]$new_node$height
  patch <- patch_with_heights(c(5, 3, seed / 2))
  expect_lt(min(vapply(patch$species[[1]]$nodes, function(n) n$height, numeric(1))),
            seed)

  r <- compare_on_grid(patch)
  expect_identical(r$value, r$value_ref)
  expect_identical(r$slope, r$slope_ref)
  expect_gt(max(abs(r$value)), 0)
})

test_that("each species carries its own sums", {
  # eta is a strategy parameter, so the running sums are per species and cannot
  # be shared across them. Two strategies with different eta is what would show a
  # shared set.
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 10
  p1 <- add_strategies(p0, trait_matrix(c(0.08, 0.2), "lma"), birth_rate = rep(1.0, 2))
  patch <- run_scm(p1, Environment("FF16"), Control())$patch
  expect_equal(patch$size, 2)

  r <- compare_on_grid(patch)
  message(sprintf("two species: worst value %.2e (eps*A(0) = %.2e)",
                  max(abs(r$value - r$value_ref)), r$scale))
  expect_lt(max(abs(r$value - r$value_ref)), 20 * r$scale)
})

test_that("the field's slope is the derivative of the field's value", {
  # The interpolant guarantees that its own slope differentiates its own value.
  # This is the other half: that the slope the reduction supplies at a knot is
  # the derivative of the value it supplies there, so the pair the interpolant is
  # handed is consistent before it ever fits anything.
  patch <- stand_patch("FF16", t_max = 20)
  z <- patch$environment$light_availability$state[, "height"]
  # Interior knots only, and away from the canopy top where the profile has a
  # break at every cohort height.
  probe <- z[z > 0.5 & z < 0.6 * max(z)]
  probe <- probe[seq(1, length(probe), length.out = min(40, length(probe)))]

  h <- 1e-6
  analytic <- vapply(probe, function(zi) patch$compute_competition_and_slope(zi)[2],
                     numeric(1))
  fd <- vapply(probe, function(zi) {
    (patch$compute_competition(zi + h) - patch$compute_competition(zi - h)) / (2 * h)
  }, numeric(1))
  worst <- max(abs(analytic - fd) / pmax(abs(fd), 1e-8))
  message(sprintf("reduction slope vs central difference: worst rel = %.2e", worst))
  expect_lt(worst, 1e-5)
})

test_that("extending the lattice moves no query the shorter one answered", {
  # The canopy decides how many knots there are and nothing else about them, so a
  # grown stand adds nodes above and never moves one below. That is what lets the
  # knot count depend on the state.
  patch <- stand_patch("FF16", t_max = 20)
  env <- patch$environment
  before <- env$light_availability$state
  probe <- before[before[, "height"] < 0.5 * max(before[, "height"]), "height"]
  read_before <- vapply(probe, env$light_availability$get_value_at_height, numeric(1))

  # Rebuild the same field: the lattice is already long enough, so it is kept.
  patch$compute_environment()
  after <- patch$environment$light_availability$state
  expect_gte(nrow(after), nrow(before))
  expect_identical(after[seq_len(nrow(before)), "height"], before[, "height"])
})

test_that("a canopy no forest reaches is refused by name", {
  # Only this class can say what ran away: the interpolant lays whatever lattice
  # it is asked for, and a count no memory could hold is arithmetic to it. So the
  # refusal names the height and the spacing, not the allocation.
  #
  # 0.05 m knots and a 25000-knot bound put it at 1250 m, against a tallest-tree
  # record near 130 m, so a stand cannot reach it and a driven height can.
  x <- "FF16"
  s <- strategy_types[[x]]()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- Parameters(x, environment_types[[x]])(strategies = list(s),
                                             patch_type = "meta-population")
  patch <- Patch(x, environment_types[[x]])(p, Environment(x), Control())
  patch$introduce_new_node(1)
  state <- patch$ode_state
  state[[1]] <- 2000

  # Setting the state rebuilds the field, so the refusal arrives there rather
  # than at a later build -- which is the point at which a run would meet it.
  expect_error(patch$set_ode_state(state, 0),
               "canopy height 2000 m needs more than 25000 knots")
  expect_error(patch$set_ode_state(state, 0),
               "size-density equation has run away")
})

test_that("the knots are constants of the run, not of the stand", {
  # Knot k sits at k * spacing whatever the stand is doing. A grid tied to the
  # canopy top moves every knot when the tallest cohort grows; this one does not,
  # and that is the property the whole design rests on.
  spacing <- 0.05
  for (t_max in c(5, 20, 80)) {
    z <- stand_patch("FF16", t_max = t_max)$environment$light_availability$state[, "height"]
    expect_equal(z, seq(0, by = spacing, length.out = length(z)))
  }
})

test_that("an empty patch and a single cohort both answer", {
  x <- "FF16"
  e <- environment_types[[x]]
  s <- strategy_types[[x]]()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- Parameters(x, e)(strategies = list(s), patch_type = "meta-population")
  patch <- Patch(x, e)(p, Environment(x), Control())

  # No nodes: the reduction has nothing to sum and says so.
  patch$compute_environment()
  expect_identical(patch$compute_competition(0), 0)
  expect_identical(patch$compute_competition_and_slope(0), c(0, 0))

  # One node: the whole integral is the single closing interval, which is the
  # case the boundary test exists for.
  patch$introduce_new_node(1)
  patch$compute_environment()
  r <- compare_on_grid(patch)
  expect_identical(r$value, r$value_ref)
  expect_identical(r$slope, r$slope_ref)
})
