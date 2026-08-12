# Rung 4: two species, two cohorts each, over a trajectory.
#
# Two species with one cohort each catch a reduction sum that collapses the
# species and not per-species accumulation; two cohorts of one species catch the
# accumulation and not the collapse. Four nodes catch both and are still small
# enough to reason about by hand.
#
# New over rung 3: the stage recursion, step-size handling, accumulation across
# cohorts and species, and record-once-sweep-many.
#
# Here the Jacobian can no longer be formed, so the contraction returns -- with
# its blind spot stated. A contraction against one direction masks a column whose
# true value is near zero, so it is paired with the zero classification of the
# floor and with coordinate directions on a declared shortlist, each chosen
# because it has a documented structural route to reading exactly zero.

ladder_shortlist <- function() {
  # Each of these reaches a census through a reduction as well as through a
  # plant, so each can read zero for a structural reason rather than an
  # ecological one.
  c("k_I", "eta", "a_l1", "a_l2", "a_st3", "recruitment_decay")
}

test_that("the structural checks hold before any reference is consulted", {
  # These need no reference, so they fail earlier and more cheaply than the
  # contraction -- and a failure here invalidates the contraction anyway.
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")
  result <- ladder_gradient_or_skip(stand)

  expect_equal(dim(result$gradient), c(3L, 88L))
  expect_true(all(is.finite(result$gradient)))

  # The trait adjoint accumulates over every cohort, every stage and every step,
  # so a trait treated as one input per cohort returns a fixed fraction of the
  # right answer with the correct sign and no error raised. Two species with two
  # cohorts each is the smallest stand where that fraction is not one.
  expect_length(census_trait_names_tf24(stand), 88L)
})

test_that("the replay reaches the run's own census", {
  # Before any derivative the reference reports means anything, it has to be a
  # derivative of the same function. The tangent run replays the recorded step
  # SIZES rather than the times -- a size differenced back out of two recorded
  # times is not the size that was taken, since fl(fl(t + h) - t) != h -- and
  # rather than a controller of its own, which would differentiate the
  # controller. Landing on the run's own census is what says the replay worked.
  #
  # It is also every trajectory check's floor: the reference reports derivatives
  # at the state it actually reached.
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")
  result <- ladder_gradient_or_skip(stand)
  ladder_report_margin("the replay reaches the run's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)

  # Non-vacuity of the seed: an unseeded replay must report no derivative, or a
  # zero tangent would be evidence of nothing about a column that reads zero.
  quiet <- ladder_trajectory_tangent(
    stand, numeric(length(colnames(result$gradient))))
  expect_true(all(quiet$tangent == 0))
})

test_that("the gradient contracts to a tangent of the same trajectory", {
  # The reference traverses both reductions, the stage recursion and the
  # introduction boundary, and none of the transposes under test are on its path.
  # Here the Jacobian can no longer be formed, so the contraction returns.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  columns <- colnames(result$gradient)

  direction <- ladder_seeds(length(columns),
                            scale = ladder_block_scale(result$gradient[1, ]))
  reference <- ladder_trajectory_tangent(stand, direction)$tangent
  observed <- as.vector(result$gradient %*% direction)
  expect_length(reference, nrow(result$gradient))

  # Normalised by the largest single contribution rather than by their sum. The
  # sum cancels -- that cancellation is the contraction's blind spot, and it is
  # why the coordinate directions below are run beside it -- so dividing by it
  # would report cancellation as error.
  term <- max(abs(sweep(result$gradient, 2, direction, `*`)))
  ladder_report_margin("the contraction against the trajectory tangent",
                       max(abs(reference - observed)) / term,
                       3 * ladder_introduction_residual())
})

test_that("each shortlisted trait's own column is refereed", {
  # A contraction against one direction masks a column whose true value is near
  # zero, so every shortlisted parameter is asked for as a coordinate direction:
  # one replay, one exact column of the census-by-trait Jacobian.
  #
  # Two bounds, because the columns fall into two classes and one bound would
  # excuse the worse of them. A column a widening carries is held to the looser
  # one, and which columns those are is declared rather than discovered from the
  # residual.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  columns <- colnames(result$gradient)
  present <- columns[ladder_bare_traits(columns) %in% ladder_shortlist()]
  expect_length(present, 10L)

  for (name in present) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    peak <- max(abs(result$gradient[, name]))
    expect_gt(peak, 0)
    borne <- ladder_bare_traits(name) %in% ladder_introduction_borne_traits()
    ladder_report_margin(
      paste("column", name, if (borne) "(introduction-borne)" else ""),
      max(abs(reference - result$gradient[, name])) / peak,
      if (borne) ladder_introduction_residual()
      else ladder_trajectory_agreement())
  }
})

test_that("no shortlisted trait reads exactly zero", {
  # A contraction against one direction masks a column whose true value is near
  # zero. Coordinate directions on the shortlist are what cover that blind spot,
  # and every entry on the list is there because it has a route to reading zero
  # for a structural reason.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  # Every species' column for every shortlisted parameter, because a column that
  # reads zero for one species only is exactly the collapse the shortlist is
  # watching for. The crown shape has no column at all and drops out here.
  columns <- colnames(result$gradient)
  present <- columns[ladder_bare_traits(columns) %in% ladder_shortlist()]
  expect_length(present, 10L)

  peak <- vapply(present, function(n) max(abs(result$gradient[, n])), numeric(1))
  message("\n  shortlisted trait magnitudes:")
  for (n in present) {
    message(sprintf("    %-20s %.6e", n, peak[[n]]))
  }
  expect_true(all(peak > 0))
})

test_that("a trait is one input read by every cohort, not one input per cohort", {
  # Because the strategy is shared, a single trait is one input. Treating each
  # cohort as a distinct input yields a gradient that is a fixed fraction of the
  # right answer with the correct sign and no error raised, and the fraction is
  # the cohort count. So the signature is a ratio, and adding a cohort at fixed
  # ecology is what exposes it.
  fewer <- ladder_stand_two_by_two()
  result_fewer <- ladder_gradient_or_skip(fewer)

  more <- local({
    p <- ladder_parameters(c("fast", "slow"))
    p$node_schedule_times <- list(c(0, 0.31, 0.63), c(0, 0.41))
    ladder_run(p)
  })
  result_more <- stand_gradient(more)

  # Refining the schedule at fixed ecology moves the answer by quadrature, not
  # by a whole cohort count. A ratio near the node-count ratio is the signature.
  ratio <- result_more$gradient["leaf_area", "1.lma"] /
    result_fewer$gradient["leaf_area", "1.lma"]
  message(sprintf("\n  lma row, four nodes over three: %.4f", ratio))
  expect_lt(abs(ratio - 1), 0.5)
})

test_that("the two species do not collapse into one scalar", {
  # A reduction that sums a per-species parameter over every cohort of every
  # species collapses the species into one number. This is undetectable with one
  # species, which is what earns the second, and the fixture's species differ by
  # whole factors in exactly the parameters that carry a per-species reduction
  # row.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  columns <- census_trait_names_tf24(stand)
  half <- length(columns) / 2L

  for (name in c("k_I", "a_l1", "a_l2")) {
    # Naming the two columns rather than searching for them: a search that comes
    # back with one match, or none, would skip the check that earns the second
    # species.
    at <- match(paste0(1:2, ".", name), columns)
    expect_false(anyNA(at), label = paste("both species' columns for", name))
    first <- result$gradient[, at[[1]]]
    second <- result$gradient[, at[[2]]]
    expect_false(isTRUE(all.equal(first, second)),
                 label = paste("the two species' rows for", name, "are equal"))
  }
  expect_equal(half, 44L)
})

test_that("record once and sweep many gives each metric its own recording", {
  # Clearing a tape returns its derivative-slot counter to zero, so an active
  # value constructed outside the sweep loop and read inside it refers, after the
  # first clear, to a slot that now belongs to something else. The first metric
  # is correct and every later one reads unrelated storage -- and a value that
  # survives a clear registers with no dependencies, so its adjoint sweeps to
  # exactly zero, which is why the failure produces exact zeros in whole column
  # families rather than noise.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  metrics <- rownames(result$gradient)

  zero_families <- vapply(metrics, function(m) all(result$gradient[m, ] == 0),
                          logical(1))
  expect_false(any(zero_families))

  # And the positional check, which CANNOT FAIL as written and is kept only
  # because the shape it watches for is real. stand_gradient(metrics = m) computes
  # every metric in C++ and subsets in R, so "swept alone" and "swept with the
  # others" are one sweep and comparing them is a tautology. Established by
  # injection: with an active value deliberately held across the tape clear, this
  # passed and only the trajectory reference saw the defect. Closing it needs the
  # sweep itself to take a metric subset.
  for (m in metrics) {
    alone <- stand_gradient(stand, metrics = m)
    expect_identical(alone$gradient[m, ], result$gradient[m, ])
  }
})
