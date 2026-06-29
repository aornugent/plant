# TF24f individual grow-to-size trait gradient (#472 scope B, the "individuals" surface).
# A single plant grown in a FIXED environment to target sizes; the gradient is the
# finite-difference d(t*)/d(theta) and d(state at t*)/d(theta) over
# grow_individual_to_size (no resident feedback, no canopy/density -- the lightest
# surface). This is the prototype + the reference an exact-AD version must reproduce;
# the tracked collar (opt_root_psi_state) is re-evolved inside each grow, so the FD
# captures its theta-response automatically.

test_that("tf24f individual grow-to-size FD gradient is well-formed and sensible", {
  s <- TF24f_Strategy()
  ind <- TF24f_Individual(s)
  env <- Environment("TF24f")
  env$set_fixed_environment(1.0, 1e4)
  sizes <- c(2, 5, 9)
  g <- plant:::tf24f_grow_individual_to_size_gradient_fd(
    ind, sizes, "height", env, traits = c("vcmax_25", "lma", "K_s"), time_max = 400)

  # Shapes + the tracked collar is one of the returned state components.
  expect_length(g$time, length(sizes))
  expect_equal(dim(g$d_time), c(length(sizes), 3L))
  expect_true("opt_root_psi_state" %in% dimnames(g$d_state)[[2]])
  expect_true(all(is.finite(g$d_time)) && all(is.finite(g$d_state)))

  # Reconstructed stopping times match grow_individual_to_size.
  ref <- grow_individual_to_size(ind, sizes, "height", env, time_max = 400)
  expect_equal(g$time, ref$time, tolerance = 1e-6)

  # Physically-signed sensitivities of the time-to-size: costlier leaves (higher lma)
  # slow growth (d t* / d lma > 0); higher stem conductance (K_s) speeds it
  # (d t* / d K_s < 0). Checked at the largest target where the signal is clear.
  expect_gt(g$d_time[length(sizes), "lma"], 0)
  expect_lt(g$d_time[length(sizes), "K_s"], 0)
})

test_that("tf24f individual grow gradient is strategy-guarded", {
  ind <- FF16_Individual(FF16_Strategy())
  env <- FF16_Environment(); env$set_fixed_environment(1.0, 1e4)
  expect_error(
    plant:::tf24f_grow_individual_to_size_gradient_fd(ind, c(2, 5), "height", env),
    "TF24f strategy only")
})
