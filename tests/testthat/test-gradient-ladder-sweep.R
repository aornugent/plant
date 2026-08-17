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

test_that("a shut point answers, and the branch it took is what decides", {
  # This block used to assert the opposite. Dried to half its moisture the
  # one-cohort patch reaches a genuine hydraulic shutdown, and the boundary
  # declined to answer for it -- so the test recorded the refusal and the fact
  # that its message named the kind.
  #
  # The shut branches answer now. What survives unchanged is the reason the old
  # check gave for reading the KIND rather than a residual: the marginal profit
  # returns a hard sentinel zero in a no-flow state, and no residual test can tell
  # that from stationarity. The classification still has to come from the exit the
  # solve took -- it is what routes this state away from an envelope step it has
  # no stationary point for.
  interior <- ladder_patch_one()
  shutdown <- ladder_patch_shutdown()

  # Both answer, and the shut one is not answered by accident: it is a different
  # derivation, reached because the leaf reports the branch.
  expect_silent(invisible(ladder_block_value_tf24(interior, 1L)))
  expect_silent(invisible(ladder_block_value_tf24(shutdown, 1L)))
})

test_that("a shut point's rows are numbers, and every one of them is finite", {
  # The old form of this check asserted that no vector was produced at all, on
  # the grounds that a refused row coming back as zeros is indistinguishable from
  # an answer. That hazard has not gone away -- it has moved: the rows here ARE
  # numbers now, so what has to hold is that they are finite and that the zeros
  # among them are the structural ones.
  #
  # Where a refusal is still reachable, `test-gradient-parity.R` carries the
  # absence-never-a-number check against a driver that genuinely refuses.
  shutdown <- ladder_patch_shutdown()
  n <- shutdown$ode_size
  seed <- ladder_seeds(n, scale = ladder_block_scale(shutdown$ode_rates))
  got <- ladder_rhs_adjoint_tf24(shutdown, seed)$state
  expect_true(is.numeric(got))
  expect_true(all(is.finite(got)))
  # Non-vacuity: a shut leaf still moves the stand, so a row of zeros would mean
  # the state had been answered for by not being read.
  expect_gt(sum(abs(got) > 0), 0)
})

test_that("both output kinds answer at a shut point, which they did not before", {
  # Report 08 asks the two output kinds to refuse independently: the profit row
  # survives every degeneracy except a jump of the argmax and an undefined
  # objective, and the uptake row is the one that ceases to exist.
  #
  # At a shut point the question is now moot in the direction that matters --
  # BOTH answer. The rows are a substituted constant plus the bound's movement,
  # and neither output kind is missing. So this measures that all four seeds come
  # back rather than that all four are refused together.
  #
  # ⚠️ The independence requirement is NOT discharged by this. The fixture that
  # would separate the two kinds is a fold, and nothing here reaches one; what
  # changed is only that this fixture no longer refuses either.
  shutdown <- ladder_patch_shutdown()
  nm <- ladder_rate_names(shutdown)
  n <- shutdown$ode_size
  kinds <- list("size only, height" = "height",
                "size only, heartwood" = "mass_heartwood",
                "water-coupled, a layer" = "environment_1",
                "water-coupled, total uptake" = "environment_9")
  for (label in names(kinds)) {
    seed <- replace(numeric(n), match(kinds[[label]], nm), 1)
    got <- ladder_rhs_adjoint_tf24(shutdown, seed)$state
    message(sprintf("  %-28s answered", label))
    expect_true(all(is.finite(got)))
  }
})
