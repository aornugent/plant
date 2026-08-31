# The sweep computes the same thing however it is decomposed.
#
# Every assertion here is an exact identity at tolerance zero, so none of them
# needs a reference and none of them is a statement about run length: composition
# over steps is associative or it is not, a sweep is repeatable or it is not, and a
# permutation of the metric order changes the numbers or it does not. That is why
# they run on the short fixture -- the same five introductions over a quarter of the
# steps -- and why nothing in this file was re-blessed to put them there.
#
# Their whole assurance is the non-vacuity clause beside each one: prove the two
# paths computed something before believing they agree.

test_that("two consecutive sweeps of one recording are bit-identical", {
  # The check that the forward replay of introductions leaves the system where
  # it found it. A run that has to be repeatable is a run whose replay is not
  # consuming the state it replays.
  # The short fixture: this asserts an exact identity, which run length cannot
  # weaken, and it is one of the three checks that drove this file's runtime.
  stand <- ladder_stand_introductions_short()
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
  # The short fixture, for the same reason: every assertion below is an exact
  # identity, so a shorter run tests the same claim at a quarter of the sweeps.
  stand <- ladder_stand_introductions_short()
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
  # into one range per width. Composition over steps is therefore associative,
  # and splitting a range must give the whole sweep BIT FOR BIT -- tolerance is
  # exactly zero, and no reference is needed because this is a property the
  # implementation either has or does not.
  #
  # What it catches is anything carried across a step boundary that is not the
  # adjoint. The trait accumulator accumulates by design, but the block
  # workspace, the tape, the knot adjoints and the strategy templates all live
  # across steps, and a split forces a clean re-entry at the cut.
  # The short fixture keeps all five introductions and a quarter of the steps, and the
  # comparison below is bit-identity, so nothing here is re-blessed.
  stand <- ladder_stand_introductions_short()
  unsplit <- census_trait_gradient_tf24(stand)
  whole <- do.call(rbind, unsplit$gradient)
  unsplit_ranges <- unsplit$ranges

  introductions <- ladder_introduction_rows(stand$store_trajectory())
  expect_gte(length(introductions), 3L)

  # A cut names the row the descent resumes at, counted from one. A introduction is
  # already such a row, and so is the row below it -- carrying the adjoint across
  # an introduction leaves the descent there -- so a cut that adds a range is a step
  # row that is neither, which is why the two either side of an introduction are two
  # rows apart rather than adjacent. The boundary case is stated at the end.
  interior <- floor((introductions[[2]] + introductions[[3]]) / 2)
  points <- list("an interior step" = interior,
                 "one step below an introduction" = introductions[[2]] - 2,
                 "one step above an introduction" = introductions[[2]] + 1,
                 "all three at once" = c(interior, introductions[[2]] - 2,
                                         introductions[[2]] + 1))
  for (name in names(points)) {
    cut <- census_trait_gradient_split_tf24(stand, points[[name]])
    split <- do.call(rbind, cut$gradient)
    ranges <- cut$ranges
    # Non-vacuity, and it is not decoration: a split landing ON a range
    # boundary is outside every range's interior and cuts nothing, so the
    # equality below would hold between two identical sweeps. The range count is
    # what says the cut happened.
    expect_gt(ranges, unsplit_ranges)
    message(sprintf("  %-28s %2.0f ranges against %2.0f", name, ranges,
                    unsplit_ranges))
    expect_identical(split, whole)
  }

  # And the boundary case stated rather than left as a trap: naming an introduction
  # itself requests a split no range contains.
  on_boundary <- census_trait_gradient_split_tf24(stand, introductions[[2]])
  expect_equal(on_boundary$ranges, unsplit_ranges)
  expect_identical(do.call(rbind, on_boundary$gradient), whole)
})
