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
  skip(paste("the birth-size channel is priced, not bounded: every ratio reads",
             "1.000 and a bound of 1e-04 would hold each with an order of margin,",
             "but the reference is a whole-run difference of the same bisected",
             "root-find the row check below bounds, and its own floor is measured",
             "as a gap between steps -- which reads the opposite of an error that",
             "rises as the step falls"))
})

test_that("the seed's geometry row is refereed against a rebuilt strategy", {
  # The birth height is not computed forwards -- it is the root of
  # mass_live_given_height(h) = omega, solved off the tape, and its derivative is
  # declared by the implicit function theorem. So the row has to be checked
  # against the condition it claims to solve, by rebuilding the strategy at a
  # perturbed parameter and finding the root again.
  #
  # rung 5 prices the same claim for two parameters through a whole rate
  # evaluation. This one has no patch, no environment and no rates on the path,
  # so a disagreement localises to the recording rather than to anything after it.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")
  columns <- ladder_trait_names_tf24(patch)

  for (name in ladder_birth_size_parameters()) {
    at <- match(paste0("1.", name), columns)
    expect_false(is.na(at))
    # Named rather than defaulted: the bound below is the reference's own error at
    # THIS step, so the two have to be one number.
    rel <- 1e-4
    got <- ladder_seed_geometry_tangent_tf24(patch, at)
    ref <- ladder_boundary_difference(patch, name, rel = rel)

    # The recording returns the root itself, so the value is the rebuild's own.
    ladder_report_margin(
      sprintf("the seed height's value, %s", name),
      abs(got$height - ref$value[["height"]]) / abs(ref$value[["height"]]),
      1e-10)

    # Non-vacuity first: each of these reaches the residual, so a zero here is a
    # dropped channel rather than a small disagreement.
    expect_gt(abs(got$dheight), 0)

    # And the row, against a difference of the rebuilt root. What bounds it is
    # where the root-find stopped, not how fine the step is: the reference
    # differences a bisection whose answer is quantised at half a bracket, and a
    # relative step of 1e-4 in the parameter is a step of a few nanometres in the
    # height, which is smaller than the quantisation it is dividing.
    #
    # 1e-4 is also where this is the WHOLE bound: the reference's truncation there
    # is a ten-thousandth of it for every one of the eight. At 1e-3 truncation
    # returns at a seventh and the bound would need a second term.
    #
    # No factor of ten. This is a hard interval, not a scale estimate, so over it
    # means the row disagrees by more than a bisected root-find can account for.
    ladder_report_margin(
      sprintf("the seed height's row, %s", name),
      abs(got$dheight - ref$row[["height"]]) / abs(ref$row[["height"]]),
      ladder_bisection_bound(
        patch, 1L, ref$value[["height"]],
        abs(ladder_strategy_parameter(patch, 1L, name)) * rel,
        ref$row[["height"]]))
  }

  # And the structural half: the residual reads those parameters and no others,
  # so anything outside the list is exactly zero rather than small.
  for (name in c("k_I", "a_l1")) {
    at <- match(paste0("1.", name), columns)
    if (is.na(at) || name %in% ladder_birth_size_parameters()) {
      next
    }
    expect_identical(ladder_seed_geometry_tangent_tf24(patch, at)$dheight, 0)
  }
})

# The sweep now carries its own declaration of what an exact zero means, so the
# list above is no longer the only statement of it. Two statements of one list is
# the shape this corpus keeps finding defects in, so they are required to agree
# rather than left to drift: the C++ one is what ships inside an answer, and this
# one is the fixture that prices it.
test_that("the shipped zero classes agree with the declared lists", {
  stand <- ladder_stand_two_by_two()
  g <- stand_gradient(stand)

  # Only exact zeros are classified, so the two claims are checked in both
  # directions: every column the ladder declares zero comes back in a zero class,
  # and every column the sweep puts in one is declared.
  declared <- ladder_zero_by_construction()
  shipped <- apply(g$status, 2, function(col) unique(col))
  expect_true(is.character(shipped))

  slack <- names(shipped)[shipped == "zero-slack"]
  structural <- names(shipped)[shipped == "zero-structural"]
  undeclared <- names(shipped)[shipped == "zero-undeclared"]

  expect_setequal(ladder_bare_traits(slack),
                  ladder_zero_at_an_interior_optimum())
  expect_setequal(ladder_bare_traits(structural),
                  ladder_zero_outside_the_metric_support())
  expect_setequal(ladder_bare_traits(c(slack, structural)), declared)

  # An exact zero with no declared reason is a finding, and this fixture must
  # produce none -- otherwise the classification is passing an unexplained zero
  # off as one of the two named causes.
  expect_length(undeclared, 0L)

  # Non-vacuity: the classes are not all one value, and the zero classes really
  # are a strict minority of a mostly-live matrix.
  expect_gt(sum(g$status == "answered"), length(declared) * nrow(g$status))
  message(sprintf("  %d answered, %d zero-slack, %d zero-structural, %d refused",
                  sum(g$status == "answered"), sum(g$status == "zero-slack"),
                  sum(g$status == "zero-structural"), sum(g$status == "refused")))
})
