# Rung 3: one cohort, and the Jacobian formed entirely.
#
# Use the strongest oracle the size permits, and the size permits an exhaustive
# one exactly once -- here, at the bottom. Nothing above this is small enough to
# form whole, so nothing above it is ever checked this well again.
#
# Do not contract. A dot-product identity yields one number per run: an error in
# one cell enters as a seed component times that cell times an input component,
# so it hides behind a small seed, it cancels against other cells, and when it
# does fail it localises to nothing. At this size the whole object can be formed,
# so forming it is strictly better -- a failure names an output row and an input
# column rather than reporting that a sum disagreed.
#
# The unit is the adjoint of one right-hand-side evaluation at one stage. No
# solver, no trajectory.

ladder_block_or_skip <- function(patch, node = 1L) {
  out <- tryCatch(list(value = ladder_block_value_tf24(patch, node)),
                  error = function(e) e)
  if (inherits(out, "error")) {
    testthat::skip(paste("the cohort block does not record at an active scalar:",
                         conditionMessage(out)))
  }
  out
}

test_that("the recorded block has the shape the design is costed on", {
  patch <- ladder_patch_one()
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  outputs <- ladder_block_output_names_tf24(patch)

  # Six own states, the field's knot values and slopes, the soil potentials and
  # the traits in; six rates, the density rate and one draw per layer out. Both
  # counts move with the configuration and neither is a property of the model;
  # the shape is.
  n_knot <- (length(inputs) - 6L - 5L - 44L) / 2L
  expect_equal(length(inputs), 6L + 2L * n_knot + 5L + 44L)
  expect_equal(length(outputs), 6L + 1L + 5L)
  expect_equal(head(inputs, 6L),
               c("height", "mortality", "fecundity", "area_heartwood",
                 "mass_heartwood", "storage"))
  expect_equal(tail(outputs, 5L), paste0("uptake_", 1:5))
  message(sprintf("\n  block shape: %d inputs (%d knots) by %d outputs",
                  length(inputs), n_knot, length(outputs)))
})

test_that("the block's forward and reverse Jacobians agree entry by entry", {
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")
  ladder_block_or_skip(patch)

  inputs <- ladder_block_input_names_tf24(patch, 1L)
  outputs <- ladder_block_output_names_tf24(patch)
  forward <- ladder_block_jacobian_forward_tf24(patch, 1L)
  reverse <- ladder_block_jacobian_reverse_tf24(patch, 1L)
  dimnames(forward) <- dimnames(reverse) <- list(outputs, inputs)

  # Non-vacuity: a Jacobian of a block that returned nothing is not evidence.
  expect_true(all(is.finite(forward)))
  expect_true(all(is.finite(reverse)))
  expect_gt(sum(forward != 0), length(outputs))

  # Each entry is measured against its own output row's scale, not against
  # itself. A cell whose true value is a hundred orders below its row carries no
  # information about the row, and a per-cell relative measure would make the
  # tolerance a statement about that cell's round-off rather than about the
  # Jacobian.
  row_scale <- pmax(apply(abs(forward), 1, max), .Machine$double.xmin)
  residual <- abs(forward - reverse) / row_scale
  worst <- which(residual == max(residual), arr.ind = TRUE)[1, ]
  message(sprintf("  worst cell: %s / %s  forward %.6e  reverse %.6e",
                  outputs[[worst[[1]]]], inputs[[worst[[2]]]],
                  forward[worst[[1]], worst[[2]]],
                  reverse[worst[[1]], worst[[2]]]))

  ladder_report_margin("block Jacobian, forward against reverse",
                       max(residual), 10 * ladder_forward_floor(patch))
})

test_that("the block's structure is what every cost argument assumes", {
  # The matrix formed above is the five-way classification, measured rather than
  # asserted. Asserting it directly costs a rank computation on a matrix that
  # already exists, and this is the one place the structure is checked rather
  # than read.
  patch <- ladder_patch_one()
  ladder_block_or_skip(patch)
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  outputs <- ladder_block_output_names_tf24(patch)
  j <- ladder_block_jacobian_forward_tf24(patch, 1L)
  dimnames(j) <- list(outputs, inputs)

  is_light <- grepl("^light_", inputs)
  is_value <- grepl("^light_value_", inputs)
  is_slope <- grepl("^light_slope_", inputs)
  is_soil <- grepl("^psi_soil_", inputs)
  uptake <- grepl("^uptake_", outputs)

  numerical_rank <- function(m, tol = 1e-9) {
    if (all(m == 0)) return(0L)
    s <- svd(m)$d
    sum(s > tol * max(s))
  }

  # The field reaches the physiology through one scalar, so the whole knot block
  # is rank one on this coordinate. On the height coordinate the density rate
  # would evaluate the physiology a second time at a displaced height and make
  # its own field query, and the block would have rank at least two.
  light_block <- j[, is_light, drop = FALSE]
  message(sprintf("  light block rank %d over %d columns",
                  numerical_rank(light_block), sum(is_light)))
  expect_equal(numerical_rank(light_block), 1L)

  # Per-layer uptake reads its own layer's potential and the collar, so the
  # explicit part is diagonal. What is left after removing the diagonal is the
  # argmax channel, and it is rank one because the collar is one number.
  soil_block <- j[uptake, is_soil, drop = FALSE]
  expect_equal(dim(soil_block), c(5L, 5L))
  message(sprintf("  uptake-by-potential block rank %d",
                  numerical_rank(soil_block)))
  off_diagonal <- soil_block
  diag(off_diagonal) <- 0
  expect_equal(numerical_rank(off_diagonal), 1L)

  # Density enters no cohort rate: it reaches the world only through the two
  # reductions. It is not among the block's inputs at all, which is the strongest
  # form of that claim.
  expect_false(any(grepl("log_density", inputs)))

  # The value and slope channels are separate inputs and both must be live. An
  # identity posed on values alone passes with the entire slope channel absent.
  ladder_expect_moves(j[, is_value], matrix(0, nrow(j), sum(is_value)),
                      "the field's value channel")
  ladder_expect_moves(j[, is_slope], matrix(0, nrow(j), sum(is_slope)),
                      "the field's slope channel")
})

test_that("every trait the block reads has a column, or is refused by name", {
  # An absent row is an exact zero that reads as an ecological finding. The block
  # is where a trait's row either exists or does not, so this is where the
  # classification is cheapest to make.
  patch <- ladder_patch_one()
  ladder_block_or_skip(patch)
  inputs <- ladder_block_input_names_tf24(patch, 1L)
  j <- ladder_block_jacobian_forward_tf24(patch, 1L)
  traits <- inputs[!grepl("^light_|^psi_soil_", inputs)][-(1:6)]
  columns <- j[, match(traits, inputs), drop = FALSE]
  peak <- apply(abs(columns), 2, max)
  names(peak) <- traits

  refused <- ladder_refused_by_name()
  declared <- ladder_zero_by_construction()
  unaccounted <- names(peak)[peak == 0 & !(names(peak) %in% c(refused, declared))]
  if (length(unaccounted) > 0) {
    message("\n  traits with no row in the block: ",
            paste(unaccounted, collapse = " "))
  }
  expect_length(unaccounted, 0L)
})

# ---- the reduction transposes, where the reference shares no code ------------
#
# The checks below run on the uncrossed fixture, and that is a concession rather
# than a choice: the crossed one refuses, and the test directly above owns that
# refusal. Once the transposes test the abscissa rather than the height, these
# should move to the crossed fixture, which is where the design says every run
# belongs.

test_that("a crossed stand transposes, because the forward model runs it", {
  # Reserve-gated growth lets a younger cohort overtake an older one, so on the
  # birth-date coordinate crossing is more common, not less. The quantity that
  # must stay monotone here is the abscissa, and the abscissa is not height -- so
  # a guard written as a test of height ordering refuses exactly the stands the
  # forward model handles correctly.
  #
  # Both failures are silent and they are silent in opposite directions: a
  # reduction that integrates a crossed grid without guarding has neighbouring
  # trapezia cancelling instead of accumulating, and a guard on the wrong
  # quantity stops on a stand that is fine.
  crossed <- ladder_patch_two_by_two(cross = TRUE)
  ladder_require_regime(crossed, "patch")
  nodes <- ladder_nodes(crossed)
  expect_gt(ladder_crossing_count(nodes), 0L)

  # The forward model runs it, which is the half of the claim that has to hold
  # before the other half is a defect rather than a preference.
  rates <- crossed$ode_rates
  expect_true(all(is.finite(rates)))

  seed <- ladder_seeds(crossed$ode_size,
                       scale = ladder_block_scale(rates))
  expect_no_error(ladder_rhs_adjoint_tf24(crossed, seed))
})

test_that("one right-hand-side evaluation transposes to its own tangent", {
  # Here the reference and the object under test are genuinely disjoint: the
  # tangent traverses the forward reductions and the adjoint traverses their
  # transposes. Every reduction defect the design fears is a wrong entry in this
  # comparison.
  patch <- ladder_patch_two_by_two(cross = FALSE)
  ladder_require_regime(patch, "patch")
  ladder_block_or_skip(patch)

  n <- patch$ode_size
  jacobian <- ladder_rhs_state_jacobian_forward_tf24(patch)
  expect_equal(dim(jacobian), c(n, n))
  expect_true(all(is.finite(jacobian)))

  # Seeds from a fixed generator, block-normalised, because soil moisture is of
  # order a third and heartwood mass is in kilograms and an unnormalised inner
  # product is a test of the largest block alone.
  rates <- patch$ode_rates
  seed <- ladder_seeds(n, scale = ladder_block_scale(rates))

  expected <- as.vector(crossprod(jacobian, seed))
  observed <- ladder_rhs_adjoint_tf24(patch, seed)$state

  scale <- max(abs(expected))
  residual <- max(abs(expected - observed)) / scale
  worst <- which.max(abs(expected - observed))
  message(sprintf("  worst state entry %d: tangent %.6e  adjoint %.6e",
                  worst, expected[[worst]], observed[[worst]]))
  ladder_report_margin("right-hand-side transpose, per state entry",
                       residual, 1e3 * ladder_forward_floor(patch))
})

test_that("the trait rows that arise inside a reduction arrive", {
  # Four parameters reach a census through the field build as well as through a
  # plant: the extinction coefficient, the crown shape, and the two allometric
  # constants that set leaf area at fixed height. A transpose written through a
  # structure holding only size and density slots has no route to a trait
  # accumulator however correct its arithmetic, and those rows then read exactly
  # zero.
  patch <- ladder_patch_two_by_two(cross = FALSE)
  ladder_block_or_skip(patch)

  n <- patch$ode_size
  rates <- patch$ode_rates
  seed <- ladder_seeds(n, scale = ladder_block_scale(rates))

  jacobian <- ladder_rhs_trait_jacobian_forward_tf24(patch)
  expected <- as.vector(crossprod(jacobian, seed))
  observed <- ladder_rhs_adjoint_tf24(patch, seed)$trait
  names(expected) <- names(observed) <- ladder_trait_names_tf24(patch)

  # Non-vacuity first: the reduction's own contribution has to be non-zero in
  # the reference, or its absence in the sweep proves nothing.
  reduction_borne <- c("k_I", "eta", "a_l1", "a_l2")
  present <- intersect(reduction_borne, names(expected))
  expect_gt(length(present), 0L)
  expect_true(all(abs(expected[present]) > 0),
              label = "the tangent's reduction-borne rows are themselves zero")

  scale <- max(abs(expected))
  residual <- max(abs(expected - observed)) / scale
  gap <- names(expected)[abs(expected - observed) > 1e-6 * scale]
  if (length(gap) > 0) {
    message("\n  trait rows that disagree: ", paste(unique(gap), collapse = " "))
  }
  ladder_report_margin("right-hand-side transpose, per trait",
                       residual, 1e3 * ladder_forward_floor(patch))
})

test_that("the recording does not grow with the stand", {
  # Peak is one cohort's rates at one stage, so it is flat in the cohort count
  # and in the number of differentiation targets: seeding more traits adds
  # registered slots, not recorded operations. Either climbing with the cohort
  # count is the tape leaking.
  one <- ladder_patch_one()
  four <- ladder_patch_two_by_two()
  ladder_block_or_skip(one)

  seed_one <- rep(1, one$ode_size)
  seed_four <- rep(1, four$ode_size)
  a <- ladder_rhs_adjoint_tf24(one, seed_one)
  b <- ladder_rhs_adjoint_tf24(four, seed_four)
  message(sprintf("\n  recording size: %.0f at 1 cohort, %.0f at 4",
                  a$block_recording_size, b$block_recording_size))
  expect_equal(b$block_recording_size, a$block_recording_size)
  expect_equal(b$block_sweeps, 4)
})
