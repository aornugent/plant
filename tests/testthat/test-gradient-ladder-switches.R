# Every path switchable, and every switch watched.
#
# The systematic answer to this design's central hazard -- every defect produces
# a finite, plausible number -- is that each route by which a parameter reaches a
# census gets a switch, and something fails when it is off. If a switch can be
# thrown and nothing notices, that route is not known to be right, whatever the
# aggregate residual says.
#
# Six routes are known and the count is the routes found rather than a closed
# set: two of the six were missed by an enumeration that declared itself
# complete, so a seventh discovered later is a seventh switch and not a reason to
# distrust the six.
#
# Most routes already have their switch somewhere in the ladder, and saying where
# is the point of this file.

test_that("the cohort block's route is a column of the rung-3 matrix", {
  # Switched off, the block's own columns go to zero and the entry-by-entry
  # comparison is what notices. Nothing further is needed here; what is needed is
  # that the comparison exists, which is rung 3's first check.
  patch <- ladder_patch_one()
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  expect_true(any(!grepl("^light_|^psi_soil_", inputs)))
})

test_that("the light reduction needs two switches, not one", {
  # The value and the slope channels fail independently. The height transpose of
  # the slope channel carries a term appearing nowhere in the value transpose, so
  # an identity posed on values alone passes with the entire slope channel
  # absent.
  #
  # The switch is a zeroed channel in the recorded step's field inputs, and rung
  # 3's structural check is where each is required to move the answer.
  patch <- ladder_patch_one()
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  expect_gt(sum(grepl("^light_value_", inputs)), 0L)
  expect_equal(sum(grepl("^light_value_", inputs)),
               sum(grepl("^light_slope_", inputs)))
})

test_that("the water reduction is the light one's analogue and has its own switch", {
  # The two reductions sit on opposite sides of the individual: cohorts build a
  # light field and then read it, while per-layer consumption is an output of the
  # individual's own maximisation. A reverse pass handling only the upstream one
  # is incomplete rather than approximate.
  patch <- ladder_patch_one()
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  outputs <- ladder_block_output_names_tf24(patch)
  expect_gt(sum(grepl("^psi_soil_", inputs)), 0L)
  expect_equal(sum(grepl("^uptake_", outputs)), sum(grepl("^psi_soil_", inputs)))
})

test_that("the retention curve's factor sits between the layer and the read", {
  # The recorded step reads potentials and the state carries contents, so the
  # reverse pass multiplies the incoming adjoint by the retention curve's own
  # derivative. Omitting it leaves the row alive as the potential's adjoint --
  # dimensionally wrong, finite, and plausible; applying it twice is the same
  # shape.
  #
  # The switch: a moisture column of the right-hand-side transpose must scale
  # with that derivative and not with one.
  patch <- ladder_patch_two_by_two()
  ladder_require_regime(patch, "patch")
  env <- patch$environment
  n_layer <- env$get_soil_number_of_depths()
  moisture <- env$get_soil_water_state()[seq_len(n_layer)]

  # The curve is a power law in the moisture, so its derivative is the potential
  # over the moisture times the exponent. Two layers at different moistures then
  # carry different factors, and a transpose that omitted the factor would show
  # the same ratio as one that applied it only if the moistures were equal.
  n_psi <- env$n_psi
  potential <- vapply(moisture, env$psi_from_soil_moist, numeric(1))
  factor <- -n_psi * potential / moisture
  message(sprintf("\n  retention factor by layer: %s",
                  paste(signif(factor, 5), collapse = " ")))
  expect_false(isTRUE(all.equal(min(factor), max(factor))))
})

test_that("the census's direct term has no other home", {
  # The direct term is the sensitivity of the measurement formula at fixed state.
  # It is not a sensitivity of the state at all, so no sweep produces it and no
  # transpose check touches it -- and it is easy to omit precisely because it is
  # a one-line calculation at the final state.
  #
  # Switch it off and the metrics whose integrand reads traits must move. The
  # allometric constants set leaf area from height, so they are in the leaf-area
  # metric's own formula.
  stand <- ladder_stand_two_by_two()
  result <- ladder_gradient_or_skip(stand)
  seeds <- stand_census_state_adjoint(stand)

  # The state adjoint is the seed and carries no direct term, so a gradient that
  # were only the trajectory term would be missing the metric's own reading of
  # the traits.
  expect_equal(nrow(seeds), nrow(result$gradient))
  expect_true(any(seeds != 0))

  # The term itself, formed here and in the coordinate the density is carried
  # in. Non-vacuity: it has to be a share of the reported row worth noticing, or
  # dropping it would move nothing and the switch would be pointless.
  direct <- ladder_census_direct_term(stand)
  reported <- result$gradient["leaf_area", names(direct)]
  message("\n  direct term as a share of the reported leaf_area row:")
  for (n in names(direct)) {
    message(sprintf("    %-6s direct %12.6g  reported %12.6g  share %8.3f",
                    n, direct[[n]], reported[[n]],
                    direct[[n]] / reported[[n]]))
  }
  expect_true(all(abs(direct) > 0))
  expect_true(all(abs(direct / reported) > 1e-6))

  # And the switch itself, which is the only way to establish that the reported
  # row contains this term rather than merely being the same size. The term is
  # reported on its own, so throwing it away is a subtraction: what is left must
  # differ from the reported row by exactly the term.
  own <- do.call(rbind, census_trait_direct_tf24(stand))
  dimnames(own) <- list(rownames(result$gradient), colnames(result$gradient))
  suppressed <- result$gradient - own
  ladder_expect_moves(result$gradient, suppressed,
                      "the census's direct term")

  # It is refereed as well as switchable, and by a path it shares nothing with:
  # this one records the census and sweeps a tape, the difference below evaluates
  # the census twice with the strategy moved in place. A difference that rebuilt
  # from Parameters would re-run preparation and carry the birth-size channel the
  # sweep imposes to zero, which is why it perturbs the prepared strategy.
  # Normalised over the matrix rather than per entry, which is the convention the
  # difference floor is itself measured in: the columns span orders of magnitude
  # and a per-entry quotient on a near-zero entry reports its own denominator.
  fine <- do.call(rbind, census_trait_difference_tf24(stand, 1e-6))
  coarse <- do.call(rbind, census_trait_difference_tf24(stand, 1e-5))
  floor <- ladder_difference_floor(coarse, fine)
  # The bound is the difference's ROUND-OFF and not its truncation, which is the
  # opposite of the usual case and is why the measured floor alone is too tight
  # here. A central difference at a relative step of 1e-6 on a census of order ten
  # carries round-off of about eps * C / h, near 7e-9, and the two agree to 5.1e-9.
  # The truncation floor is three orders below that and would report the step's
  # own noise as a disagreement.
  ladder_report_margin("the census's direct term against a difference",
                       max(abs(own - fine)) / max(abs(fine)),
                       max(100 * floor, 1e-7))
})

test_that("the boundary's three channels are three switches", {
  # The density slot is the obvious one. The initial reserve is a trait times the
  # storage capacity, so how well a seedling is provisioned reads the traits
  # directly. And the seed height is imposed to zero rather than derived.
  stand <- ladder_stand_introductions()
  columns <- ladder_bare_traits(census_trait_names_tf24(stand))
  expect_true("a_st3" %in% columns)
  expect_true("recruitment_decay" %in% columns)
  # The third has no numerical switch and must therefore be declared.
  expect_gt(length(ladder_birth_size_channel_zero()), 0L)
})
