# AD-9 (#14): the TF24 leaf edge inside the SCM -- PROTO-2, the highest-risk item.
# TF24 solves a per-cohort leaf hydraulic optimum whose per-area carbon profit drives
# growth. A leaf trait (vcmax_25) reaches an emergent density-dependent stand metric
# ONLY through that optimum -- it is absent from the allometric mass / area / fecundity
# model -- so the whole invasion gradient w.r.t. vcmax_25 IS the density->optimum
# cross-term: the trait moves the optimum -> growth -> cohort density and size -> the
# metric. The leaf's own solve stays a plant-local, off-tape double computation; at each
# rate evaluation its trait and size sensitivities are re-attached to the one reverse
# tape (a native linearisation about the operating point), so the sweep traverses
# density->optimum->trait rather than as a frozen harvest (which would zero the
# cross-term). The reverse-mode gradient must match a re-optimising central finite
# difference of the same frozen-canopy mutant run.

test_that("TF24 invasion gradient matches re-optimising central FD (density->optimum cross-term)", {
  skip_on_cran()  # a TF24 leaf-hydraulic SCM run + finite-difference oracle is expensive

  mpl <- 6
  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- mpl
  e <- Environment("TF24")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE
  # Tighten the leaf optimiser so the re-optimising finite difference is clean (the AD
  # side reuses the same off-tape solve for its size-feedback sensitivity).
  ctrl$GSS_tol_abs <- 1e-9
  ctrl$ci_abs_tol <- 1e-10

  pr1 <- add_strategies(p0, trait_matrix(96, "vcmax_25"), birth_rate = 20)
  # Mature early so the single invading species reproduces within the short run, giving
  # a substantial survival-weighted offspring signal to differentiate.
  pr1$strategies[[1]]$pars$hmat <- 5
  scm <- run_scm(pr1, e, ctrl)

  # AD-7's native census kernels are exercised for TF24: the LAI / biomass / basal-area
  # invasion gradients must be finite and non-zero (each carries the trait through the
  # leaf optimum into the cohort size distribution).
  cen <- invasion_gradient(scm, c("biomass", "LAI", "basal_area"), "vcmax_25")
  expect_true(all(is.finite(cen$gradient)))
  expect_true(all(abs(cen$gradient) > 0))

  # Density-dependent offspring gradient: the smooth survival-weighted stand integral.
  # vcmax_25 reaches it only through the leaf optimum, so this gradient is entirely the
  # density->optimum cross-term.
  ad <- offspring_production_gradient(scm, "vcmax_25")
  off0 <- scm$offspring_production[1]
  expect_true(is.finite(ad))

  # Cross-term is substantial: a frozen-optimum harvest (leaf profit held constant w.r.t.
  # the trait) would be identically zero, so a large gradient is the cross-term itself.
  expect_gt(abs(ad / off0), 1e-2)

  # Re-optimising central FD (Richardson-extrapolated to strip the O(h^2) bias): perturb
  # the raw vcmax_25 the AD seeds, re-run the frozen-canopy mutant, read offspring. The
  # optimum re-solves each perturbation, so a non-zero FD is the empirical proof the
  # census responds through the leaf optimum -- not a frozen harvest.
  fd_off <- function(delta) {
    p <- pr1
    p$strategies[[1]]$pars$vcmax_25 <- p$strategies[[1]]$pars$vcmax_25 + delta
    scm$run_mutant(p)
    scm$offspring_production[1]
  }
  x0 <- pr1$strategies[[1]]$pars$vcmax_25
  central <- function(h) (fd_off(h * x0) - fd_off(-h * x0)) / (2 * h * x0)
  richardson <- function(h) (4 * central(h / 2) - central(h)) / 3
  fd <- richardson(1e-4)

  expect_gt(abs(fd / off0), 1e-2)              # the optimum genuinely responds
  expect_equal(ad, fd, tolerance = 1e-4)       # AD captures the density->optimum cross-term
})

test_that("TF24f trait gradients fail loudly (coupled feedback out of v1)", {
  skip_on_cran()

  p0 <- scm_base_parameters("TF24f", "TF24_Env")
  p0$max_patch_lifetime <- 2
  e <- Environment("TF24")
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE
  pr1 <- add_strategies(p0, trait_matrix(96, "vcmax_25"), birth_rate = 20)
  scm <- run_scm(pr1, e, ctrl)

  # TF24f's tracked-leaf coupling is stiff; the fixed-step replay drifts, so the entry
  # gates with a clear v1-scope error rather than returning a wrong number.
  expect_error(invasion_gradient(scm, "biomass", "vcmax_25"), "TF24f")
  expect_error(stand_gradient(scm, "biomass", "vcmax_25"), "TF24f")
})
