# TF24f stand census reconstruction -- the R0 GATE (#472 scope B, build-order step 1).
# Per notes/tf24-stand-gradient-scope.md, TF24f is the right target for the stand CENSUS
# gradients (its tracked-collar leaf evaluation is analytic / IFT-able, so the census
# number density needs no curvature harvest). Before any reverse-mode tape (R1) is built,
# this gate proves the prerequisite §7 flags: that re-evolving each cohort's {5 demog,
# tracked collar, log_density} over the SCM's frozen Cash-Karp schedule -- driving the
# real TF24f leaf at the tracked collar in double precision -- reproduces the SCM's stored
# stand (heights / opt_root_psi_state / log_density) and the census LAI =
# compute_competition(0). The recon is value-only (no AD), so it runs in plain R / CI.
#
# Patch lifetime 4 keeps the tallest cohort ~7 m, inside the regime where the replay is
# bit-faithful (heights / collar to ~1e-6, log_density to ~1e-6). For taller cohorts the
# growth-rate-gradient finite difference (g' = backward FD of height_dt, abs step 1e-6) is
# the fidelity-limiting term: as growth saturates with height, g' is a near-cancelling
# difference of two large crown-centre leaf solves, and at lifetime 5 (tallest ~8.4 m) the
# oldest cohort's log_density drifts ~6e-3 (LAI ~0.4%); at lifetime 8 the tracked collar
# drifts out of the leaf transport spline domain entirely (an acclimation
# schedule-sensitivity to pin before pushing to long horizons).

tf24f_small_scm <- function(H = 4L, n = 9L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
}

# Read the SCM's stored per-cohort final states. node_at() is 1-based; the ODE state
# layout is [height, mortality, fecundity, area_heartwood, mass_heartwood,
# opt_root_psi_state, offspring, log_density].
tf24f_scm_stand <- function(scm, species = 1L) {
  sp <- scm$patch$species[[species]]
  nc <- length(sp$node_times)
  list(
    heights = vapply(seq_len(nc), function(i) sp$node_at(i)$ode_state[1], numeric(1)),
    collar  = vapply(seq_len(nc), function(i) sp$node_at(i)$ode_state[6], numeric(1)),
    log_density = sp$log_densities,
    LAI = scm$patch$compute_competition(0))
}

test_that("tf24f_census_recon reproduces the SCM's collar-state stand (R0 gate)", {
  scm <- tf24f_small_scm()
  ref <- tf24f_scm_stand(scm)
  rec <- plant:::tf24f_census_recon(scm, metrics = c("LAI", "biomass", "size_moment"))

  # Per-cohort heights and the tracked collar (opt_root_psi_state) -- the replay
  # re-evolves the 6-state ODE, so matching these confirms the collar-state trajectory.
  expect_equal(rec$heights, ref$heights, tolerance = 1e-4)
  expect_equal(rec$collar,  ref$collar,  tolerance = 1e-4)
  # The census number density (the genuinely new state). g' (the SCM's backward-FD
  # growth-rate gradient) is the fidelity-limiting term, hence the looser absolute floor.
  expect_equal(rec$log_density, ref$log_density, tolerance = 1e-3)
})

test_that("tf24f_census_recon LAI reduction matches compute_competition(0)", {
  scm <- tf24f_small_scm()
  ref <- tf24f_scm_stand(scm)
  rec <- plant:::tf24f_census_recon(scm, metrics = "LAI")
  # The height-trapezium of density * k_I * area_leaf (+ pending-seed tail) is exactly
  # the SCM's competition integral at ground level.
  expect_equal(rec$values[["LAI"]], ref$LAI, tolerance = 1e-3)
  # biomass / size_moment are the same trapezium with a different per-cohort kernel;
  # they must be finite and positive on a live stand.
  rec2 <- plant:::tf24f_census_recon(scm, metrics = c("biomass", "size_moment"))
  expect_true(all(is.finite(rec2$values)) && all(rec2$values > 0))
})

test_that("tf24f_census_gradient_fd: finite frozen census Jacobian (R1 gate)", {
  # The reverse-mode census trait gradient (R1) must reproduce a two-pass FD over the
  # SAME frozen-env replay; this builds that FD reference and checks it is well-formed
  # and self-consistent (a usable prototype gradient pending the AD tape).
  scm <- tf24f_small_scm()
  tr <- c("vcmax_25", "lma", "a_l1", "K_s")
  g <- plant:::tf24f_census_gradient_fd(scm, metrics = c("LAI", "biomass", "size_moment"),
                                        traits = tr)
  expect_equal(dim(g$jacobian), c(3L, length(tr)))
  expect_true(all(is.finite(g$jacobian)))
  # The reconstructed metric values match the recon (and so the SCM census).
  expect_equal(g$values[["LAI"]], scm$patch$compute_competition(0), tolerance = 1e-3)
  # Self-consistency: a column recomputed at a 10x coarser FD step agrees to ~1% -- the
  # precision of any FD over the recon, capped by the leaf optimiser/root-find noise
  # floor (~5e-7) amplified by the step. This is the tolerance band the reverse-mode AD
  # tape will be validated within (the AD removes this FD-step sensitivity).
  g2 <- plant:::tf24f_census_gradient_fd(scm, metrics = "LAI", traits = "lma",
                                         rel_step = 1e-4)
  expect_equal(unname(g$jacobian["LAI", "lma"]), unname(g2$jacobian["LAI", "lma"]),
               tolerance = 2e-2)
})

test_that("tf24f resident census FD gradient is the total gradient (flips sign)", {
  # The resident (coupled) gradient is the TOTAL stand-level d(metric)/d(theta): the
  # canopy feedback routinely dominates and flips the sign relative to the frozen
  # (rare-mutant / invasion) gradient. Ground truth here is FD over the full SCM.
  H <- 4L
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = 9L))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  tr <- c("lma", "K_s")
  gr <- plant:::tf24f_resident_census_gradient_fd(p, Environment("TF24f"), ctlc,
          metrics = c("LAI", "size_moment"), traits = tr)
  expect_equal(dim(gr$jacobian), c(2L, length(tr)))
  expect_true(all(is.finite(gr$jacobian)))

  scm <- run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
  expect_equal(unname(gr$values[["LAI"]]), scm$patch$compute_competition(0),
               tolerance = 1e-6)
  # The defining property: feedback flips d(LAI)/d(lma) -- positive for a rare mutant
  # against the frozen canopy, negative as the whole stand's LMA rises.
  gf <- plant:::tf24f_census_gradient_fd(scm, metrics = "LAI", traits = "lma")
  expect_gt(gf$jacobian["LAI", "lma"], 0)
  expect_lt(gr$jacobian["LAI", "lma"], 0)
})

test_that("tf24f census recon is strategy-guarded", {
  # A TF24 (non-f) resident has no tracked-collar state; the harvest must reject it.
  p <- scm_base_parameters("TF24"); p$max_patch_lifetime <- 4L
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 4L, length.out = 9L))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  scm <- run_scm(p, Environment("TF24"), ctlc, refine_schedule = FALSE)
  expect_error(plant:::tf24f_census_recon(scm), "TF24f strategy only")
})
