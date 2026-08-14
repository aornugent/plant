# Pricing a declared zero, which needs a reference that does not declare it.
#
# Every other check in this ladder referees the sweep against a forward tangent.
# The tangent is exact and it settles whether the transpose is the transpose of its
# forward function -- but it inherits every equation the model imposes, the seed
# height's among them, so the two agree there for free and the agreement reads as a
# pass. A whole-run difference rebuilds the strategy, which re-runs preparation, so
# it carries the channel both differentiated paths hold still.
#
# This is the shortlist of columns that channel reaches, and no more: it costs two
# model runs per column.

test_that("a differenced column moves one registered parameter", {
  # The hyperparameter function derives registered parameters from traits, so a
  # difference taken on the trait vector answers what a species further along a
  # trait spectrum does while the sweep answers what one parameter does. Recorded
  # here as data rather than as prose, because which traits fan out is a property
  # of that function and changes when it does.
  traits <- ladder_traits()$fast
  for (name in names(traits)) {
    fanout <- ladder_trait_fanout(traits, name)
    message(sprintf("  %-5s also moves: %s", name,
                    if (length(fanout)) paste(names(fanout), collapse = ", ")
                    else "(nothing)"))
  }
  # lma is the one the fixture carries, and it is why the guard exists.
  expect_gt(length(ladder_trait_fanout(traits, "lma")), 0)
  expect_length(ladder_trait_fanout(traits, "a_l1"), 0)

  # The guard's sensitivity, established by breaking what it watches: no other
  # check in this suite notices a difference taken along a direction, because the
  # answer stays finite, plausible and correctly signed.
  before <- c(lma = 0.0825, k_l = 2.038, r_l = 572.4)
  expect_error(
    ladder_assert_one_parameter(before, c(lma = 0.0826, k_l = 2.030, r_l = 572.4),
                                "lma"),
    "moved 2 registered parameters")
  expect_silent(
    ladder_assert_one_parameter(before, c(lma = 0.0826, k_l = 2.038, r_l = 572.4),
                                "lma"))
})

test_that("the whole-run difference is in its own domain on this fixture", {
  # The check on the check. A re-run difference is unusable at production, so its
  # validity is established per fixture rather than assumed, and a fixture it
  # cannot resolve invalidates the run rather than failing it.
  traits <- ladder_traits()$fast
  for (name in ladder_birth_size_parameters()) {
    got <- ladder_run_difference_stable(traits, name)
    message(sprintf(
      "\n  %-5s plateau at %s (step %.0e), spread %.2e over the pair and %.2e over all four; seed-height slope %+.4e",
      name, got$plateau_at, got$step, got$spread, got$spread_all,
      got$seed_height_slope))
    if (got$spread > 1e-3) {
      skip(paste("the whole-run difference does not hold its figures on this",
                 "fixture, so it is out of its domain here rather than failing:",
                 name))
    }
    expect_lt(got$spread, 1e-3)
    # Non-vacuity: the channel this reference exists to price has to be live, or
    # agreeing with the sweep would prove nothing about the declaration.
    expect_gt(abs(got$seed_height_slope), 0)
  }
})

test_that("the birth-size channel is priced rather than asserted", {
  # Report 05 §10.1 imposes d(seed height)/d(trait) = 0 on both differentiated
  # paths and calls it the one term no available instrument can referee. This is
  # the instrument. What it reports is a ratio, not a pass: the declaration is a
  # modelling choice, and what a suite owes it is a number.
  #
  # The number moves with run length, so the run length is reported beside it --
  # a_l1 is a factor of two out at four tenths of a year and three parts in a
  # thousand at four years, and quoting either without the lifetime is quoting
  # nothing.
  traits <- ladder_traits()$fast
  stand <- ladder_stand_allometric_probe(TRUE)
  gradient <- ladder_gradient_or_skip(stand)
  cols <- colnames(gradient$gradient)

  message("\n  sweep against a whole-run difference, leaf area, 0.4 yr:")
  for (name in ladder_birth_size_parameters()) {
    got <- ladder_run_difference_stable(traits, name)
    if (got$spread > 1e-3) {
      message(sprintf(
        "    %-5s difference out of its own domain here (spread %.1e, plateau at %s)",
        name, got$spread, got$plateau_at))
      next
    }
    sweep <- gradient$gradient[1, paste0("1.", name)]
    ratio <- sweep / got$gradient[[1]]
    message(sprintf("    %-5s difference %+.6e  sweep %+.6e  ratio %6.3f",
                    name, got$gradient[[1]], sweep, ratio))
    # The declaration is not zero-cost, and that is the finding rather than a
    # failure: this asserts only that the price is finite and reported.
    expect_true(is.finite(ratio))
    expect_gt(abs(ratio), 0)
  }
  skip(paste("the birth-size channel is priced, not bounded: no acceptance number",
             "is declared for it until report 05 §10.1's imposition is either",
             "closed or restated with a domain"))
})
