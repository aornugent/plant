
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]


  test_that("Empty environment", {
    p <- Parameters(x, e)()
    env <- Environment(x)

    ## At this point, we should have full canopy openness, partly because
    ## the spline is just not constructed.
    expect_equal(env$get_environment_at_height(0), 1.0)
    expect_equal(env$get_environment_at_height(100), 1.0)

    ## Before the first build the field is the open one: three knots carrying the
    ## same value and no slope, which reads as 1 everywhere it is defined.
    state <- env$light_availability$state
    expect_equal(nrow(state), 3)
    expect_equal(state[, "height"], c(0, 0.5, 1))
    expect_equal(state[, "light_availability"], rep(1, 3))
    expect_equal(state[, "slope"], rep(0, 3))
  })

  test_that("Manually set environment", {
    env <- Environment(x)
    ## Now, set the light environment. The field carries a value AND a slope at
    ## every knot, so what is handed in is the heights, then the values, then the
    ## slopes -- the columns of $state laid end to end.
    hh <- seq(0, 10, length.out=101)
    light_env <- function(z) {
      exp(z/(max(hh)*2)) - 1 + (1 - (exp(.5) - 1))/2
    }
    light_env_slope <- function(z) exp(z/(max(hh)*2)) / (max(hh)*2)
    ee <- light_env(hh)
    mm <- light_env_slope(hh)

    ## The field is a value member, so reading it gives a copy: build the field
    ## on that copy and assign it back, as the setter is the way in.
    la <- env$light_availability
    la$init_interpolators(c(hh, ee, mm))
    env$light_availability <- la

    state <- env$light_availability$state
    expect_identical(state[, "height"], hh)
    expect_identical(state[, "light_availability"], ee)
    expect_identical(state[, "slope"], mm)

    ## A knot is reproduced rather than approximated, and between the knots the
    ## field is the cubic those two ends pin.
    reads <- sapply(hh, env$light_availability$get_value_at_height)
    expect_identical(reads, ee)
    hmid <- (hh[-1] + hh[-length(hh)])/2
    expect_equal(sapply(hmid, env$light_availability$get_value_at_height),
                 light_env(hmid), tolerance = 1e-8)
  })
}

test_that("resource spline value is floored at zero (#253)", {
  # A cubic between two knots can leave the interval its ends span, and a
  # resource availability that does so is non-physical. The span here starts low
  # and steeply falling and ends flat, which is what carries it below zero --
  # a value-and-slope field needs the slope to do it, where a value-fitted one
  # could do it from the values alone.
  env <- Environment("K93")

  hh <- c(0, 1, 2)
  ee <- c(1, 0.02, 0.02)
  mm <- c(-5, -5, 0)
  la <- env$light_availability
  la$init_interpolators(c(hh, ee, mm))
  env$light_availability <- la

  grid <- seq(0, 2, by = 0.01)
  floored <- vapply(grid, env$light_availability$get_value_at_height, numeric(1))

  # The scenario genuinely undershoots, otherwise the test proves nothing. Read
  # the unclamped cubic off its own coefficients rather than off the accessor,
  # which is the thing under test.
  t <- (grid[grid > 1 & grid < 2] - 1)
  a <- ee[[2]]; b <- ee[[3]]; sa <- mm[[2]]; sb <- mm[[3]]
  c2 <- 3 * (b - a) - 2 * sa - sb
  c3 <- 2 * (a - b) + sa + sb
  raw <- a + t * (sa + t * (c2 + t * c3))
  expect_true(any(raw < 0))

  # ... but the accessor never returns a negative value ...
  expect_true(all(floored >= 0))
  # ... and is a no-op wherever the field is non-negative.
  expect_identical(floored[grid <= 1], vapply(grid[grid <= 1],
    env$light_availability$get_value_at_height, numeric(1)))
  expect_true(all(floored[grid > 1 & grid < 2] == pmax(0, raw)))
})
