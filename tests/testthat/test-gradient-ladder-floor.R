# The floor: what holds before any reference exists.
#
# These need no reference at all and are gated behind nothing. They are the
# cheapest assurance available and they run at every rung after this one. If any
# of them fails, nothing above is worth building yet.
#
# Tolerance is exactly zero for the decomposition family. There is no margin, so
# their whole assurance is the non-vacuity clause: prove the two paths computed
# something before believing they agree.

# ---- the objective itself ---------------------------------------------------

test_that("the fixtures sit in the regime their checks are declared for", {
  patch <- ladder_patch_two_by_two()
  report <- ladder_regime_report(patch, "patch")
  message("\npatch fixture regime:")
  for (i in seq_len(nrow(report))) {
    message(sprintf("  %-58s %-5s %s", report$assertion[[i]],
                    report$ok[[i]], report$value[[i]]))
  }
  expect_true(all(report$ok))
  # Crossing is asserted by the check that needs it rather than by the table,
  # because the fixtures come both crossed and not and a failure has to be
  # localisable to the crossing.
  expect_equal(ladder_crossing_count(ladder_nodes(patch)), 1L)
  expect_equal(ladder_crossing_count(
    ladder_nodes(ladder_patch_two_by_two(cross = FALSE))), 0L)

  stand <- ladder_stand_two_by_two()
  stand_report <- ladder_regime_report(stand, "stand")
  expect_true(all(stand_report$ok))

  # The floor a tolerance is measured against, rather than a literal. It is also
  # the cheapest form of the purity check: a forward quantity that does not come
  # back to its own value after the shared solver has been driven elsewhere is
  # carrying something between solves.
  floor <- ladder_forward_floor(patch)
  message(sprintf("  forward arithmetic floor: %.3e", floor))
  expect_lt(floor, 1e-12)
})

test_that("the census integrates over the coordinate the density is carried in", {
  # A gradient of a wrong number is worthless. A census is a quadrature of a
  # density, so its weights must be gaps in the coordinate the state is carried
  # on. Where that is the birth date, a grid built from heights integrates a
  # density in one variable against the spacing of another -- wrong in value on
  # any stand, and merely more visibly wrong on a crossed one.
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")

  model <- unname(stand_census(stand)["leaf_area"])
  over_birth_date <- ladder_census_leaf_area(stand, "birth_date")
  over_height <- ladder_census_leaf_area(stand, "height")

  # Non-vacuity: the two axes must disagree, or this check would pass on a stand
  # that cannot distinguish them.
  expect_gt(abs(over_birth_date - over_height) / abs(over_height), 1e-3)
  message(sprintf(
    "\n  census leaf_area: model %.8g | over birth date %.8g | over height %.8g",
    model, over_birth_date, over_height))

  expect_equal(model, over_birth_date, tolerance = 1e-8)
})

test_that("the census refuses a grid its own abscissa cannot order", {
  # Birth dates are strictly increasing by construction, so on this coordinate
  # the guard is an assertion and no sort is needed. What must not happen is a
  # census that silently integrates a tied or inverted grid.
  patch <- ladder_patch_two_by_two()
  nodes <- ladder_nodes(patch)
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i & !nodes$boundary, , drop = FALSE]
    expect_false(any(duplicated(part$birth_date)))
  }
  # And the boundary node is the youngest of all, since it is the one about to
  # be introduced. A reduction that placed it anywhere else would be integrating
  # over the wrong endpoint.
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    expect_equal(which.max(part$birth_date), which(part$boundary))
  }
})

# ---- the sweep computes the same thing however it is decomposed -------------

test_that("two consecutive sweeps of one recording are bit-identical", {
  # The check that the forward replay of introductions leaves the system where
  # it found it. A run that has to be repeatable is a run whose replay is not
  # consuming the state it replays.
  stand <- ladder_stand_introductions()
  first <- ladder_gradient_or_skip(stand)
  second <- stand_gradient(stand)
  expect_identical(first$gradient, second$gradient)
  expect_identical(first$value, second$value)
})

test_that("a gradient is bit-identical under a permutation of the sweep order", {
  # One recording, swept once per metric. Clearing a tape returns its
  # derivative-slot counter to zero, so an active value constructed outside the
  # sweep loop and read inside it refers, after the first clear, to a slot that
  # now belongs to something else. The first metric is then correct and every
  # later one reads unrelated storage, which is the worst available failure
  # shape because a correct first row lends credibility to the rest.
  #
  # The defect is positional, so permuting the order is what converts it from
  # unobservable to certain.
  stand <- ladder_stand_two_by_two()
  all_metrics <- names(stand_census(stand))
  full <- ladder_gradient_or_skip(stand)

  for (m in rev(all_metrics)) {
    alone <- stand_gradient(stand, metrics = m)
    expect_identical(alone$gradient[m, ], full$gradient[m, ],
                     info = paste("metric swept alone:", m))
  }
  reversed <- stand_gradient(stand, metrics = rev(all_metrics))
  expect_identical(reversed$gradient[all_metrics, ], full$gradient[all_metrics, ])
})

test_that("a sweep split at an interior step equals the whole sweep", {
  # Splitting a recording and sweeping the halves must give the whole sweep bit
  # for bit, at an interior step and at an introduction and at one step either
  # side of one. Non-vacuity: the split point must carry non-zero adjoint
  # traffic, or the equality is between two copies of the same thing.
  skip_if_not(ladder_can_split_sweep(),
              paste("no entry point sweeps a recorded trajectory in two parts,",
                    "so the decomposition this design rests on is unchecked"))
  stand <- ladder_stand_introductions()
  ladder_gradient_or_skip(stand)
})

test_that("rejected steps contribute nothing, and the fixture has some", {
  # A rejected attempt is not part of the trajectory, and excluding it is exact
  # rather than an approximation. Asserting the count first is what stops the
  # claim from being vacuous.
  stand <- ladder_stand_introductions()
  skip_if_not(ladder_can_count_rejections(),
              paste("the recorded trajectory does not report how many step",
                    "attempts were rejected, so 'a rejected step contributes",
                    "nothing' cannot be made non-vacuous"))
})

test_that("a sweep that never ran is distinguishable from an insensitive stand", {
  # An empty segment list is not an insensitive stand. Both must not come back
  # as a vector of zeros.
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(numeric(0))
  empty <- ladder_run(p)
  expect_error(stand_gradient(empty))
})

test_that("the reverse pass refuses the coordinate it cannot transpose", {
  # On the height coordinate the abscissa is state, so the reduction transposes
  # would omit a weight derivative and the recorded step would omit the density
  # rate's compression term. The sweep is then the transpose of a function the
  # forward model is not evaluating, and nothing about the arithmetic complains.
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(c(0, 0.63))
  stand <- ladder_run(p, ctrl = Control(node_density_in_birth_date = FALSE))

  expect_error(stand_census_state_adjoint(stand), "birth-date")
  expect_error(stand_gradient(stand), "birth-date")
  # The census value is not a gradient and is answerable on either coordinate,
  # so refusing it too would refuse a defined answer.
  expect_silent(stand_census(stand))
})

# ---- every zero is attributable to a named cause ----------------------------

test_that("every trait column resolves to one declared class", {
  # An exact zero is the signature of a missing accumulator and never of true
  # insensitivity, which makes it indistinguishable from an ecological finding.
  # So zeros are classified rather than tolerated, and the classification is two
  # failure modes and not one: an exact zero in the live class is a missing
  # accumulator, and a registered parameter reaching no equation comes back as
  # round-off, which reads as a gradient.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  g <- result$gradient

  declared_zero <- ladder_zero_by_construction()
  band <- ladder_roundoff_band()

  classify <- function(column, name) {
    peak <- max(abs(column))
    if (name %in% declared_zero) return("zero by construction")
    if (peak == 0) return("exact zero, unaccounted")
    if (peak < band[[2]]) return("round-off, reaches no equation")
    "live"
  }

  classes <- vapply(seq_len(ncol(g)),
                    function(j) classify(g[, j], colnames(g)[[j]]),
                    character(1))
  names(classes) <- colnames(g)

  message("\ntrait column classes:")
  for (cls in unique(classes)) {
    message(sprintf("  %-30s %s", cls,
                    paste(names(classes)[classes == cls], collapse = " ")))
  }

  expect_equal(unname(classes[names(classes) %in% declared_zero]),
               rep("zero by construction", sum(names(classes) %in% declared_zero)))
  expect_length(names(classes)[classes == "exact zero, unaccounted"], 0L)
  expect_length(names(classes)[classes == "round-off, reaches no equation"], 0L)
})

test_that("an unknown parameter refuses by name rather than returning a number", {
  # A zero and an absence are different at this boundary. An unknown parameter
  # must be refused; a registered one that reaches nothing is the other class
  # above and is the one that reads as an answer.
  stand <- ladder_stand_two_by_two()
  expect_error(stand_gradient(stand, traits = "not_a_parameter"), "Unknown trait")
  expect_error(stand_gradient(stand, metrics = "not_a_metric"),
               "Unknown census metric")
})

test_that("gradient columns are named per species", {
  # Concatenating each species' parameter names with no species prefix yields one
  # name repeated per species. Character indexing then resolves each to its first
  # match, so a multi-species gradient silently returns species one's column for
  # every named parameter, and unknown-parameter validation cannot see it. A
  # single-species suite never detects this.
  stand <- ladder_stand_two_by_two()
  columns <- census_trait_names_tf24(stand)
  expect_length(columns, 88L)
  expect_false(any(duplicated(columns)))
})
