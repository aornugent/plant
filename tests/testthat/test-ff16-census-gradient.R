# AD-7 (#12): the FF16 invasion census gradients (LAI / biomass / basal area) against
# a re-optimising central finite difference. Each metric is the descending-height
# trapezium of the active cohorts (heights and number density both carry the trait
# derivative) reducing the strategy's own allocation model; the reverse-mode Jacobian
# must match a central FD of the same frozen-canopy mutant run. Exercises the compiled
# invasion_gradient path in plant.so.

test_that("FF16 invasion census gradients match re-optimising central FD", {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  e <- Environment("FF16")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE

  pr1 <- add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = 1)
  scm <- run_scm(pr1, e, ctrl)

  metrics <- c("LAI", "biomass", "basal_area")
  traits  <- c("lma", "hmat")

  ad <- invasion_gradient(scm, metrics, traits)
  expect_true(all(is.finite(ad$gradient)))
  expect_equal(dim(ad$gradient), c(length(metrics), length(traits)))

  # Re-optimising central FD of the same census metric: perturb the raw strategy
  # parameter the AD seeds, re-run the frozen-canopy mutant and read the metric value
  # (the identical kernel the AD differentiates). Richardson-extrapolate to strip the
  # O(h^2) truncation bias so the comparison isolates the AD's exactness.
  fd_val <- function(param, delta) {
    p <- pr1
    p$strategies[[1]]$pars[[param]] <- p$strategies[[1]]$pars[[param]] + delta
    scm$run_mutant(p)
    invasion_gradient(scm, metrics, param)$value
  }
  central <- function(param, x0, h) (fd_val(param, h * x0) - fd_val(param, -h * x0)) / (2 * h * x0)
  richardson <- function(param, x0, h) (4 * central(param, x0, h / 2) - central(param, x0, h)) / 3

  for (tr in traits) {
    x0 <- pr1$strategies[[1]]$pars[[tr]]
    fd <- richardson(tr, x0, 1e-4)
    for (mt in metrics) {
      expect_equal(unname(ad$gradient[mt, tr]), unname(fd[[mt]]), tolerance = 1e-4)
    }
  }
})
