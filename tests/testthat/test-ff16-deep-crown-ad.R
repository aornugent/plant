# The FF16 invasion census gradients (LAI / biomass / basal area) through the
# compiled plant.so path. Each metric is the descending-height trapezium of the
# active mutant cohorts reducing the strategy's own allocation model -- the crown
# integral #540 approximated with a frozen weighted sum now runs as the real
# moving-node quadrature over [0, height]. Correctness is a re-optimising central FD
# of the same frozen-canopy run; regression is the recorded baseline. AD needs the
# installed DLL, so this defers to R CMD check and skips under load_all.

test_that("FF16 invasion census gradients run on the compiled path and match FD", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  pr <- gradient_fixture_parameters()
  scm <- gradient_fixture_scm()
  traits <- gradient_fixture_traits
  metrics <- gradient_fixture_metrics$resident  # the three census metrics

  ad <- invasion_gradient(scm, metrics, traits)
  expect_equal(dim(ad$gradient), c(length(metrics), length(traits)))
  expect_true(all(is.finite(ad$gradient)))

  # Re-optimising central FD of the same census metric on the frozen-canopy mutant:
  # perturb the raw strategy parameter the AD seeds, re-run, read the metric value
  # (the identical kernel the AD differentiates), Richardson-extrapolated.
  fd_val <- function(param, delta) {
    p <- pr
    p$strategies[[1]]$pars[[param]] <- p$strategies[[1]]$pars[[param]] + delta
    scm$run_mutant(p)
    invasion_gradient(scm, metrics, param)$value
  }
  central <- function(param, x0, h) (fd_val(param, h * x0) - fd_val(param, -h * x0)) / (2 * h * x0)
  richardson <- function(param, x0, h) (4 * central(param, x0, h / 2) - central(param, x0, h)) / 3

  for (tr in traits) {
    x0 <- pr$strategies[[1]]$pars[[tr]]
    fd <- richardson(tr, x0, 1e-4)
    for (mt in metrics) {
      expect_equal(unname(ad$gradient[mt, tr]), unname(fd[[mt]]), tolerance = 1e-4)
    }
  }

  base <- read_gradient_baseline()
  expect_matches_gradient_baseline(
    ad$gradient, base$jacobians$invasion$gradient[metrics, traits],
    base$fingerprint, "invasion census gradient")
})
