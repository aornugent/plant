
test_that("Gradients agree", {
  ## Test the simple finite differencing gradient function.
  gradient_fd_forward <- function(f, x, dx) {
    (f(x + dx) - f(x)) / dx
  }
  gradient_fd_centre <- function(f, x, dx) {
    (f(x + dx/2) - f(x - dx/2)) / dx
  }
  gradient_fd_backward <- function(f, x, dx) {
    (f(x - dx) - f(x)) / (-dx)
  }

  f <- function(x) x*x - 3*x + 1

  dx <- 0.001
  x <- 1

  ## Computing f(x)
  expect_identical(test_gradient_fd1(f, x, dx, 1), gradient_fd_forward(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, 0), gradient_fd_centre(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, -1), gradient_fd_backward(f, x, dx))

  ## Providing f(x)
  expect_identical(test_gradient_fd1(f, x, dx, 1, f(x)), gradient_fd_forward(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, 0, f(x)), gradient_fd_centre(f, x, dx))
  expect_identical(test_gradient_fd1(f, x, dx, -1, f(x)), gradient_fd_backward(f, x, dx))


  d <- 1e-6
  r <- 4L
  method_args <- list(d=d, eps=d)
  expect_equal(test_gradient_richardson(f, x, d, r), numDeriv::grad(f, x, method.args=method_args))
})

test_that("the reverse pass refuses a coordinate it cannot transpose", {
  # On the height coordinate the abscissa is state, so the reduction transposes
  # would omit a weight derivative and the recorded step would omit the density
  # rate's compression term. The sweep is then the transpose of a function the
  # forward model is not evaluating, and nothing about the arithmetic complains.
  # Both reverse-mode entry points must refuse, and refuse before the trait
  # gradient re-runs the model to record its trajectory.
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                      hyperpar = TF24_hyperpar, birth_rate = list(20))
  ctrl <- Control(node_density_in_birth_date = FALSE)
  scm <- run_scm(p, Environment("TF24"), ctrl, collect = FALSE)

  expect_error(stand_census_state_adjoint(scm), "birth-date")
  expect_error(stand_gradient(scm), "birth-date")

  # The census value itself is not a gradient and is answerable on either
  # coordinate, so refusing it too would refuse a defined answer.
  expect_silent(stand_census(scm))
})
