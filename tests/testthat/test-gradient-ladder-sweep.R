# Whether the sweep runs at all, and on which channels.
#
# This is the ladder's one gate. Everything that needs a trajectory sweep skips
# with the message this test fails with, so a sweep that cannot run is one red
# line naming its cause rather than the same cause repeated down every file.

test_that("the sweep runs on a stand in the declared regime", {
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")
  blocked <- ladder_sweep_blocked(stand)
  if (!is.null(blocked)) {
    message("\n  the sweep is blocked: ", blocked)
  }
  expect_null(blocked)
})

test_that("the two output kinds of the leaf refuse independently", {
  # The profit row survives every degeneracy except a jump of the argmax and an
  # undefined objective; the uptake row is the one that ceases to exist. So a
  # metric seeded only on size states survives what kills a water-coupled one,
  # and refusing them as a pair throws away the surviving metric for nothing.
  #
  # A boundary that refuses the whole sweep because one output kind has no rows
  # has made exactly that trade.
  stand <- ladder_stand_two_by_two()
  blocked <- ladder_sweep_blocked(stand)
  skip_if(is.null(blocked), "the sweep is not blocked, so there is nothing to scope")
  expect_false(grepl("per-layer uptake", blocked, fixed = TRUE),
               label = paste("the water channel's absence blocks every metric,",
                             "including ones it cannot reach:", blocked))
})

test_that("refusal is metric-level and not per-parameter", {
  # A sum has no defined value with an undefined term, so one refused operating
  # point anywhere in one metric's sweep makes that metric's whole gradient
  # undefined -- not the cohort's column, not the parameter's entry. Metrics are
  # independent of each other, so refusal is metric-level and no wider.
  #
  # A partly-populated gradient vector is the shape that reads as an answer.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  for (m in rownames(result$gradient)) {
    row <- result$gradient[m, ]
    expect_false(any(is.na(row)) && !all(is.na(row)),
                 label = paste("metric refused for some parameters only:", m))
  }
})
