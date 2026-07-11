# The FF16 invasion (offspring) gradient through the compiled plant.so path: seed
# a trait, run the real FF16_Strategy<active>/SCM, read offspring_production_gradient.
# This replaces the out-of-tree sourceCpp kernel replica (#540's FF16ProdPars
# composition, now deleted) with the production class the model actually runs.
# Correctness is pinned by a re-optimising central FD of the same frozen-canopy run;
# regression by the recorded gradient-baseline.rds oracle. AD needs the installed
# DLL, so this defers to R CMD check and skips under load_all.

test_that("FF16 invasion offspring gradient runs on the compiled path and matches FD", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  pr <- gradient_fixture_parameters()
  scm <- gradient_fixture_scm()
  traits <- gradient_fixture_traits

  ad <- offspring_production_gradient(scm, traits)
  expect_named(ad, traits)
  expect_true(all(is.finite(ad)))

  # Re-optimising central FD of the same frozen-canopy mutant: perturb the raw
  # strategy parameter the AD seeds, re-run run_mutant, read offspring. Richardson
  # extrapolation strips the O(h^2) bias so the check isolates the AD's exactness.
  fd_offspring <- function(param, delta) {
    p <- pr
    p$strategies[[1]]$pars[[param]] <- p$strategies[[1]]$pars[[param]] + delta
    scm$run_mutant(p)
    scm$offspring_production[1]
  }
  central <- function(param, x0, h) {
    (fd_offspring(param, h * x0) - fd_offspring(param, -h * x0)) / (2 * h * x0)
  }
  richardson <- function(param, x0, h) {
    (4 * central(param, x0, h / 2) - central(param, x0, h)) / 3
  }
  for (param in traits) {
    x0 <- pr$strategies[[1]]$pars[[param]]
    expect_equal(unname(ad[param]), richardson(param, x0, 1e-4), tolerance = 1e-4)
  }

  base <- read_gradient_baseline()
  expect_matches_gradient_baseline(
    ad, base$jacobians$invasion$gradient["offspring_production", traits],
    base$fingerprint, "invasion offspring gradient")
})

test_that("the gradient baseline round-trips and catches an injected perturbation", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  base <- read_gradient_baseline()
  fresh <- flatten_gradient_jacobians(gradient_fixture_jacobians())
  recorded <- flatten_gradient_jacobians(base$jacobians)

  # Round-trip: a fresh run reproduces the recorded oracle.
  expect_true(gradient_values_agree(fresh, recorded, base$fingerprint)$ok)

  # Injected perturbation: bump one recorded value well past the noise floor and
  # confirm the comparison flags it, so the oracle genuinely gates.
  bad <- recorded
  bad[1] <- bad[1] * 1.01 + 1e-6
  expect_false(gradient_values_agree(fresh, bad, base$fingerprint)$ok)
})
