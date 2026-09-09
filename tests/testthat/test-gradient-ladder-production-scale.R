# What a production-length stand violates, recorded rather than enforced.
#
# Two of the regime's assertions cannot hold at the length a user runs: heights
# pile up at the canopy so they stop being non-commensurate, and reserves
# saturate so the growth gate's slope collapses. Both are properties of the
# model at that length and not of the fixture, so enforcing them skips the run --
# and a suite that skips where it is meant to scale reports green for not having
# run. The record is the deliverable here: what a scale fixture violates is what
# the next regime has to cover.
#
# Off unless PLANT_LADDER_SCALE names a lifetime. One stand costs minutes to
# build, so this belongs in a job rather than in the loop.

test_that("a scale fixture records its regime instead of skipping", {
  stand <- ladder_stand_scale()
  nodes <- ladder_nodes(stand)

  # Must not skip: that is the whole property being added here.
  report <- ladder_require_regime(stand, "scale")
  expect_false(any(report$enforced))

  message(sprintf("\n  scale fixture: %d nodes at %g years",
                  length(nodes$height), stand$parameters$max_patch_lifetime))
  for (i in seq_len(nrow(report))) {
    message(sprintf("  %-52s %-5s %s", report$assertion[[i]],
                    if (report$ok[[i]]) "ok" else "NO", report$value[[i]]))
  }
  violated <- report$assertion[!report$ok]
  message("  violated at this length: ",
          if (length(violated)) paste(violated, collapse = "; ") else "none")

  # Non-vacuity. A scale fixture that came back the size of a localisation one
  # would record a regime nobody is scaling into.
  expect_gt(length(nodes$height), 20)
})

test_that("the two fixture classes divide the regime between them", {
  # The localisation fixtures carry the adversariality a scale fixture cannot:
  # non-commensurate heights are what makes a swapped scatter target visible, and
  # they are only available where the heights are written rather than reached.
  small <- ladder_regime_report(ladder_patch_one(), "patch")
  commensurate <- "no two heights equal or in a small-integer ratio"
  expect_true(small$ok[small$assertion == commensurate])

  stand <- ladder_stand_scale()
  large <- ladder_regime_report(stand, "scale")
  expect_false(large$ok[large$assertion == commensurate])
})
