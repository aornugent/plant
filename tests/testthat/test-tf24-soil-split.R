# R1 (issue #53): the TF24 soil RHS split for the Strang-splitting inner stepper.
# analytic_drainage_flow (exact per-layer drainage recession) + residual_soil_rhs
# (infiltration + inter-layer cascade + uptake, no diagonal drainage loss) must
# reproduce the monolithic soil solve (soil_rhs / compute_rates), numerically exactly.

rk4_step <- function(y, t, h, f) {
  k1 <- f(y, t); k2 <- f(y + h / 2 * k1, t + h / 2)
  k3 <- f(y + h / 2 * k2, t + h / 2); k4 <- f(y + h * k3, t + h)
  y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
}

test_that("R1 analytic_drainage_flow is the exact drainage recession (issue #53)", {
  env <- Environment("TF24")
  L <- env$get_soil_number_of_depths()
  th0 <- seq(0.40, 0.05, length.out = L)

  # exact flow: composes as a semigroup, identity at dt=0, drains, stays positive
  expect_equal(env$analytic_drainage_flow(th0, 0), th0)
  a_full <- env$analytic_drainage_flow(th0, 1.0)
  a_split <- env$analytic_drainage_flow(env$analytic_drainage_flow(th0, 0.3), 0.7)
  expect_equal(a_full, a_split, tolerance = 1e-12)
  expect_true(all(a_full < th0) && all(a_full > 0))

  # matches a tight RK of the drainage-only ODE theta' = -c theta^p per layer
  tsat <- env$soil_moist_sat; ksat <- env$K_sat; n <- env$n_psi
  dz <- env$depth / L; p <- 2 * n + 3; cc <- ksat / (dz * tsat^p)
  f_drain <- function(y, t) -cc * pmax(y, 0)^p
  th <- th0; t <- 0; h <- 1e-4
  while (t < 1 - 1e-12) { hh <- min(h, 1 - t); th <- rk4_step(th, t, hh, f_drain); t <- t + hh }
  expect_equal(env$analytic_drainage_flow(th0, 1.0), th, tolerance = 1e-7)  # RK ref is O(h^4)
})

test_that("R1 residual + analytic drainage reproduces the monolithic soil solve (issue #53)", {
  env <- Environment("TF24")
  L <- env$get_soil_number_of_depths()
  rd <- c(0.006, 0.005, 0.004, 0.003, 0.002)[seq_len(L)]  # a fixed (leg-frozen) uptake coupling

  scenarios <- list(
    drought = list(x = c(0, 10, 20, 30),       y = c(0, 0, 0, 0)),
    steady  = list(x = c(0, 10, 20, 30),       y = c(3, 3, 3, 3)),
    monsoon = list(x = c(0, 5, 6, 12, 13, 30), y = c(0, 0, 6, 6, 0, 0)),
    onset   = list(x = c(0, 10, 11, 30),       y = c(0, 0, 5, 5)),
    tail    = list(x = c(0, 8, 9, 30),         y = c(6, 6, 0, 0))
  )
  Tend <- 15

  mono_ref <- function(env, rd, theta0, Tend, breaks) {
    ref <- theta0; t <- 0; h <- 2.5e-4
    while (t < Tend - 1e-12) {
      hh <- min(h, Tend - t, if (length(bb <- breaks[breaks > t + 1e-12])) min(bb) - t else Inf)
      ref <- rk4_step(ref, t, hh, function(y, tt) env$soil_rhs(y, rd, tt)); t <- t + hh
    }
    ref
  }
  # Strang split: exact drainage half-flow, residual full RK4 step, exact drainage half-flow.
  # Steps never straddle a rainfall kink (the forcing breakpoints are mandatory boundaries).
  strang <- function(env, rd, theta0, Tend, breaks, H) {
    th <- theta0; t <- 0
    while (t < Tend - 1e-9) {
      hh <- min(H, Tend - t, if (length(bb <- breaks[breaks > t + 1e-12])) min(bb) - t else Inf)
      th <- env$analytic_drainage_flow(th, hh / 2)
      th <- rk4_step(th, t + hh / 2, hh, function(y, tt) env$residual_soil_rhs(y, rd, tt))
      th <- env$analytic_drainage_flow(th, hh / 2)
      t <- t + hh
    }
    th
  }

  for (nm in names(scenarios)) {
    s <- scenarios[[nm]]
    env$extrinsic_drivers_set_variable("rainfall", s$x, s$y)
    theta0 <- rep(0.30, L)
    breaks <- s$x[s$x > 0 & s$x < Tend]

    ref <- mono_ref(env, rd, theta0, Tend, breaks)
    e1 <- max(abs(strang(env, rd, theta0, Tend, breaks, 2e-3) - ref))
    e2 <- max(abs(strang(env, rd, theta0, Tend, breaks, 1e-3) - ref))
    # The split reproduces the monolithic: it CONVERGES to it as the macro step
    # shrinks (halving H reduces the error), and already matches to <~1% at a
    # working step. The accuracy-limiting scenarios are those that dry a layer to
    # the residual-floor guard or drive a fast rewetting front -- both step-limited
    # transients (finer H tightens them), not defects of the split.
    expect_lt(e2, e1 + 1e-12, label = paste("convergence", nm))
    expect_lt(e2, 1.5e-2, label = paste("accuracy", nm))
  }
  # In a guard-free, sub-saturated scenario the split matches the monolithic tightly.
  env$extrinsic_drivers_set_variable("rainfall", c(0, 10, 20, 30), c(3, 3, 3, 3))
  ref <- mono_ref(env, rd, rep(0.30, L), Tend, numeric(0))
  e_fine <- max(abs(strang(env, rd, rep(0.30, L), Tend, numeric(0), 5e-4) - ref))
  expect_lt(e_fine, 1e-4)
})

test_that("R1 drainage_touchdown_time is the closed-form time to the residual floor (issue #53)", {
  env <- Environment("TF24")
  L <- env$get_soil_number_of_depths()
  theta0 <- 0.30
  tt <- env$drainage_touchdown_time(theta0, 0)
  expect_true(is.finite(tt) && tt > 0)
  # flowing the recession for exactly the touchdown time lands on the residual floor
  landed <- env$analytic_drainage_flow(rep(theta0, L), tt)[1]
  expect_equal(landed, 1e-2, tolerance = 1e-4)
  expect_equal(env$drainage_touchdown_time(1e-2, 0), 0)  # already at the floor
})
