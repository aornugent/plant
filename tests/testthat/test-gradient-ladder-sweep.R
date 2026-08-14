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

test_that("a refusal is reachable, and it names the branch it took", {
  # The three tests above skip on a healthy stand, so until there is a fixture that
  # actually refuses they assert the shape of an answer rather than the shape of a
  # refusal. This is that fixture, and it needs no injected fault: dried to half its
  # moisture the one-cohort patch reaches a genuine hydraulic shutdown, which is one
  # of the five kinds of operating point and one the boundary declines to answer for.
  interior <- ladder_patch_one()
  shutdown <- ladder_patch_shutdown()

  # Non-vacuity first, and it is the whole of this check's validity: the same call
  # on the same fixture at full moisture has to succeed, or "it refused" would be
  # about the instrument rather than about the state.
  expect_silent(invisible(ladder_block_value_tf24(interior, 1L)))

  refusal <- tryCatch({ ladder_block_value_tf24(shutdown, 1L); NA_character_ },
                      error = function(e) conditionMessage(e))
  expect_false(is.na(refusal))
  message(sprintf("\n  the refusal: %s", refusal))

  # Reported by the branch taken, never inferred from a residual. The marginal
  # profit returns a hard sentinel zero in a no-flow state and no residual test can
  # tell that from stationarity, so the classification has to come from the exit the
  # solve took -- and it does: the message names the kind.
  expect_match(refusal, "interior optimum")
  expect_match(refusal, "shutdown")
})

test_that("a refusal is an absence and never a number", {
  # An exact zero is the signature of a missing accumulator in this design and never
  # of true insensitivity, so a refused row coming back as zeros would be the worst
  # available outcome: indistinguishable from an answer. What must happen instead is
  # that no vector is produced at all.
  shutdown <- ladder_patch_shutdown()
  n <- shutdown$ode_size
  seed <- ladder_seeds(n, scale = ladder_block_scale(shutdown$ode_rates))
  got <- tryCatch(ladder_rhs_adjoint_tf24(shutdown, seed)$state,
                  error = function(e) e)
  expect_s3_class(got, "error")
  expect_false(is.numeric(got))
})

test_that("both output kinds refuse together, which is the declared scope", {
  # Report 08 asks the two output kinds to refuse independently: the profit row
  # survives every degeneracy except a jump of the argmax and an undefined
  # objective, the uptake row is the one that ceases to exist, and a metric seeded
  # only on size states should survive what kills a water-coupled one.
  #
  # They do not, and the reason is structural rather than an oversight. The refusal
  # is raised while the block is being RECORDED, which happens before any output
  # adjoint is applied, so no seed can reach a decision that has already been taken.
  # This measures that rather than asserting the requirement, because the state it
  # is measured at is a shutdown -- where both rows are a substituted constant and
  # refusing both is the documented scope, not a lost row. The fixture that would
  # separate them is a fold, and nothing here reaches one.
  #
  # If the two are ever separated, this check fails and the domain statement that
  # calls refusal metric-level needs rewriting with it.
  shutdown <- ladder_patch_shutdown()
  nm <- ladder_rate_names(shutdown)
  n <- shutdown$ode_size
  kinds <- list("size only, height" = "height",
                "size only, heartwood" = "mass_heartwood",
                "water-coupled, a layer" = "environment_1",
                "water-coupled, total uptake" = "environment_9")
  for (label in names(kinds)) {
    seed <- replace(numeric(n), match(kinds[[label]], nm), 1)
    got <- tryCatch({ ladder_rhs_adjoint_tf24(shutdown, seed); "answered" },
                    error = function(e) "refused")
    message(sprintf("  %-30s %s", label, got))
    expect_identical(got, "refused")
  }
})
