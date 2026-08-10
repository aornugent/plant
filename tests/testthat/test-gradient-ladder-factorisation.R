# The one seam, and the check that closes it.
#
# Everywhere else in the ladder the reference and the object under test are
# disjoint. Here they are not: the tangent consumes the same supplied leaf rows
# the sweep does, and that is licensed only by those rows having been refereed
# independently against the individual's own algebra.
#
# The licence is weakest exactly where it matters most. The water rows are built
# on a factorisation of the mixed second derivative of profit into two scalars
# times closed-form vectors, and the natural check on it provably cannot referee
# it: the potential-family vectors are numerically collinear, so a compensating
# pair of coefficients fits every row equally well and an error in the second is
# absorbed into the first at a ratio of about a hundred thousand. Since the
# uniform-drying direction is a near-symmetry whose true response is amplified
# fifteen- to twenty-six-fold, a one-percent error in that coefficient is a
# fifteen- to twenty-six-fold error in the quantity the ecology cares about.
#
# So a wrong coefficient sits in the sweep and in the reference together, and
# every rung above passes. This is a common-mode failure on the channel the
# design names as its entire difficulty: the carbon half is free and the water
# half is the problem.
#
# It costs one state, no gradient run and no stand, and it gates the meaning of
# everything above it -- which is why it is the first thing to close and why it
# has a file of its own.

test_that("the marginal profit's factorisation is measured, not assumed", {
  # PASS is two things:
  #
  #   1. the second coefficient checked against its own closed form, at the
  #      first coefficient held fixed. Checking the pair jointly cannot work:
  #      a joint residual is flat along the compensating direction.
  #   2. the factorisation's residual over all of the state directions -- the
  #      soil potentials, the per-layer root masses and leaf area -- at states
  #      that include a pinned operating point and a layer near the
  #      equal-potentials branch.
  #
  # Everything the check needs is a forward quantity and each is reachable: the
  # marginal profit at a held collar potential, total uptake at a held collar
  # potential, and the physiology re-set at a perturbed soil potential, root
  # network or leaf area. The residual is then three central differences per
  # direction and a two-column least squares.
  #
  # WHAT IS NOT SETTLED, and it is why this is written as a specification rather
  # than run: the factorisation's second term differentiates the sensitivity of
  # total uptake to "root mass", and whether that is the per-layer root mass in
  # the direction being differentiated or one scalar summarising the network
  # decides both the arity of the term and what the coefficient means. Writing
  # the check against the wrong reading produces a residual that fails for a
  # reason that is not the model's.
  skip(paste(
    "the factorisation of the marginal profit is unmeasured, and it is the",
    "load-bearing claim for the whole water channel: if it does not hold, the",
    "argmax channel is not two scalars times closed-form vectors and the cost",
    "argument for the design changes. Closing it needs the reading of the",
    "root-mass direction settled first."))
})

test_that("the leaf's supplied rows are refereed against its own algebra", {
  # A finite difference of the recorded step cannot referee a supplied row. The
  # grafted expression is the value plus a sum of partials times brackets that
  # are each exactly zero, so the block's forward value does not depend on a
  # grafted input at all: differencing the block returns identically zero on
  # exactly the columns a supplied row occupies, whether the row is right,
  # wrong, or absent.
  #
  # So the rows have to be checked against the solver's own algebra, or against
  # a transpose identity that needs no reference gradient -- never against a
  # difference of the step that consumes them.
  patch <- ladder_patch_one()
  inputs <- ladder_block_input_names_tf24(patch, 1L)

  # The columns a supplied row occupies, so the claim above is concrete rather
  # than general: the radiation the leaf reads and each layer's potential.
  supplied <- grepl("^psi_soil_", inputs)
  expect_gt(sum(supplied), 0L)

  skip(paste(
    "nothing compares the leaf's supplied rows against the leaf's own",
    "algebra at one solved operating point, so the rows every rung above",
    "consumes are unrefereed"))
})
