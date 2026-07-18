# R-D (#57): the TF24 soil ODE state is log-depletion zeta = ln(theta - theta_res).
# The R interface stays in theta (set/get convert), and positivity (theta > theta_res)
# is structural -- no clamp, and theta <= 0 is not representable.

test_that("R-D: TF24 soil interface stays in theta over a log-depletion state (#57)", {
  env <- Environment("TF24")
  th <- c(0.35, 0.30, 0.25, 0.18, 0.12)
  env$set_soil_water_state(th)
  expect_equal(env$get_soil_water_state(), th)              # theta in -> theta out

  # positivity is structural: sub-residual input is floored at the residual moisture,
  # and every recovered theta is strictly above theta_res (= 0.01).
  env$set_soil_water_state(rep(0.005, 5))                   # below the residual floor
  got <- env$get_soil_water_state()
  expect_true(all(got >= 0.01))
  expect_true(all(got < 0.0100001))

  # a fresh/resized env defaults to a physical theta (half-saturation), not zeta = 0.
  env$set_soil_number_of_depths(3)
  expect_equal(env$get_soil_water_state(), rep(0.428 / 2, 3))
})
