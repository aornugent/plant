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

test_that("the reference computes the model's own rates", {
  # The cheapest check in the corpus and the one everything above it rests on: a
  # reference whose VALUE disagrees with the model is a reference to a different
  # function, and every Jacobian taken from it is that function's Jacobian --
  # internally consistent, plausible, and refereeing nothing.
  #
  # It needs no derivative and no seed, so it is a floor item rather than part of
  # the harness, and it is not subsumed by any agreement between a tangent and a
  # sweep: both are built on the same rebound twin, so a twin that differs from
  # the model in value can have them agree with each other exactly.
  for (fixture in list(ladder_patch_one(),
                       ladder_patch_two_by_two(cross = FALSE),
                       ladder_patch_two_by_two(cross = TRUE))) {
    model <- fixture$ode_rates
    reference <- ladder_rhs_value_forward_tf24(fixture)
    scale <- pmax(abs(model), abs(reference))
    relative <- ifelse(scale > 0, abs(model - reference) / scale, 0)
    worst <- which.max(relative)
    names(relative) <- ladder_rate_names(fixture)
    ladder_report_margin(
      sprintf("reference value, worst at %s", names(relative)[[worst]]),
      max(relative), 1e3 * ladder_forward_floor(fixture))
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
  # The reverse pass is a backward linear recursion over recorded steps, chopped
  # into one segment per width. Composition over steps is therefore associative,
  # and splitting a segment must give the whole sweep BIT FOR BIT -- tolerance is
  # exactly zero, and no reference is needed because this is a property the
  # implementation either has or does not.
  #
  # What it catches is anything carried across a step boundary that is not the
  # adjoint. The trait accumulator accumulates by design, but the block
  # workspace, the tape, the knot adjoints and the strategy templates all live
  # across steps, and a split forces a clean re-entry at the cut.
  stand <- ladder_stand_introductions()
  ladder_gradient_or_skip(stand)
  whole <- do.call(rbind, census_trait_gradient_tf24(stand))
  unsplit_ranges <- census_adjoint_segments_tf24(stand)

  widths <- vapply(stand$store_trajectory(), function(s) length(s$state),
                   numeric(1))
  widening <- which(diff(widths) != 0)
  expect_gte(length(widening), 3L)

  # An interior step, one either side of a widening, and all three at once. The
  # step counted here is the one the split names, from one.
  interior <- floor((widening[[2]] + widening[[3]]) / 2)
  points <- list("an interior step" = interior,
                 "one step before a widening" = widening[[2]] - 1,
                 "one step after a widening" = widening[[2]] + 1,
                 "all three at once" = c(interior, widening[[2]] - 1,
                                         widening[[2]] + 1))
  for (name in names(points)) {
    split <- do.call(rbind,
                     census_trait_gradient_split_tf24(stand, points[[name]]))
    ranges <- census_adjoint_segments_tf24(stand)
    # Non-vacuity, and it is not decoration: a split landing ON a segment
    # boundary is outside every segment's interior and cuts nothing, so the
    # equality below would hold between two identical sweeps. The range count is
    # what says the cut happened.
    expect_gt(ranges, unsplit_ranges)
    message(sprintf("  %-28s %2.0f ranges against %2.0f", name, ranges,
                    unsplit_ranges))
    expect_identical(split, whole)
  }

  # And the boundary case stated rather than left as a trap: naming a widening
  # itself requests a split no segment contains.
  boundary_split <- do.call(rbind,
                            census_trait_gradient_split_tf24(stand,
                                                             widening[[2]]))
  expect_equal(census_adjoint_segments_tf24(stand), unsplit_ranges)
  expect_identical(boundary_split, whole)
})

test_that("a rejected step attempt is not a question about the gradient", {
  # This is where the ladder asked for the rejection count, and the requirement is
  # withdrawn rather than skipped. The reverse pass never sees a rejected attempt:
  # the adaptive pass resolves the schedule, one state and one step size are
  # recorded per ACCEPTED step, and the sweep runs over those. There is no
  # exclusion left to perform, so "the gradient is unchanged when rejections are
  # excluded" compares a thing with itself.
  #
  # The seam a rejected attempt does leave is in the FORWARD replay, not the
  # sweep: it moves patch state that is not ODE state -- the first-same-as-last
  # derivative carry, and anything cached on the patch -- so a pinned replay is
  # not bit-identical to the adaptive run it replays. That is what the trajectory
  # reference's own floor measures, and it is asserted where the reference is
  # used rather than here.
  stand <- ladder_stand_introductions()
  result <- ladder_gradient_or_skip(stand)
  ladder_report_margin("the replay reaches this stand's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)
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
  bare <- ladder_bare_traits(colnames(g))

  classify <- function(column, name) {
    peak <- max(abs(column))
    declared <- name %in% declared_zero
    if (declared && peak == 0) return("zero by construction")
    # A declared zero that comes back non-zero is as much a finding as an
    # undeclared zero: the declaration names a reason, and a number means the
    # reason has stopped holding.
    if (declared) return("declared zero, now live")
    if (peak == 0) return("exact zero, unaccounted")
    if (peak < band[[2]]) return("round-off, reaches no equation")
    "live"
  }

  classes <- vapply(seq_along(bare),
                    function(j) classify(g[, j], bare[[j]]), character(1))
  names(classes) <- colnames(g)

  message("\ntrait column classes:")
  for (cls in unique(classes)) {
    message(sprintf("  %-30s %s", cls,
                    paste(names(classes)[classes == cls], collapse = " ")))
  }

  expect_length(names(classes)[classes == "declared zero, now live"], 0L)
  expect_length(names(classes)[classes == "exact zero, unaccounted"], 0L)
  expect_length(names(classes)[classes == "round-off, reaches no equation"], 0L)
  # Non-vacuity: a declaration that named nothing would make the class above
  # unfalsifiable, and every column being live would make this check pass with
  # the classification absent.
  expect_gt(sum(classes == "zero by construction"), 0L)
  expect_gt(sum(classes == "live"), 0L)
})

test_that("each declared zero is zero for the cause it is declared for", {
  # Naming a cause is only worth more than tolerating a zero if the cause is
  # measured. Each cause is a claim about which rows of the right-hand side the
  # parameter moves, and the trait Jacobian of one rate evaluation answers it
  # directly: a column of zeros for a parameter claimed to move one rate, or a
  # moved rate for one claimed to move none, is a wrong declaration.
  # One species, because the cause is a property of the parameter and not of
  # which species carries it, and the trait Jacobian costs one evaluation per
  # column: forty-four here against eighty-eight on the four-node fixture, for
  # the same claim. The per-species question belongs to rungs 3 and 4.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")
  j <- ladder_rhs_trait_jacobian_forward_tf24(patch)
  columns <- ladder_trait_names_tf24(patch)
  bare <- ladder_bare_traits(columns)
  rates <- ladder_rate_names(patch)

  message("\ndeclared zeros, and the rates they move:")
  for (name in ladder_zero_by_construction()) {
    # Both species' columns, so a cause that holds for one species and not the
    # other is a failure rather than a pass on the first match.
    at <- which(bare == name)
    expect_length(at, length(patch$species))
    moved <- rates[apply(abs(j[, at, drop = FALSE]) > 0, 1, any)]
    message(sprintf("  %-16s %-46s %s", name, ladder_zero_cause(name),
                    paste(unique(moved), collapse = " ")))
    if (name %in% ladder_zero_at_an_interior_optimum()) {
      # The dry bound of a feasible interval the operating point is inside. It
      # moves nothing at all, and that is complementary slackness rather than a
      # missing row -- the same parameter carries the whole row at a pin.
      expect_length(moved, 0L)
    } else {
      # The two reproductive accumulators and nothing else. A third rate would
      # make the census's silence about the column wrong, and neither of these
      # two is read by any metric or by any equation -- which is why the same
      # column would be live on a fitness functional.
      expect_setequal(unique(moved), ladder_reproductive_rates())
    }
  }
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
