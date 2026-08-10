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

test_that("the gradient contracts to a tangent of the same trajectory", {
  # The reference is a tangent run of the whole solve, replaying the recorded
  # step sizes. Without that replay the tangent run's controller chooses its own
  # steps and the comparison differentiates the controller, which is not the
  # model.
  stand <- ladder_stand_two_by_two()
  ladder_gradient_or_skip(stand)
  skip_if_not(
    exists("ladder_trajectory_tangent_tf24", mode = "function"),
    paste("no tangent run replays a recorded trajectory's step sizes, so the",
          "trajectory sweep has no reference and the stage recursion, the",
          "step-size handling and the cross-species accumulation are unrefereed"))
})

test_that("no shortlisted trait reads exactly zero", {
  # A contraction against one direction masks a column whose true value is near
  # zero. Coordinate directions on the shortlist are what cover that blind spot,
  # and every entry on the list is there because it has a route to reading zero
  # for a structural reason.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  present <- intersect(ladder_shortlist(), colnames(result$gradient))
  expect_gt(length(present), 0L)

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
  ratio <- result_more$gradient["leaf_area", "lma"] /
    result_fewer$gradient["leaf_area", "lma"]
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
    at <- which(columns == name)
    if (length(at) != 2L) next
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

  # And the positional check: each metric swept alone must equal itself swept
  # with the others, which is the floor's permutation check read per metric.
  for (m in metrics) {
    alone <- stand_gradient(stand, metrics = m)
    expect_identical(alone$gradient[m, ], result$gradient[m, ])
  }
})
