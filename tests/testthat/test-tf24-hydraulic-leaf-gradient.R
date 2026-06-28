# TF24 HYDRAULIC leaf-trait gradients (#472 scope B, Phase F1-full): exact
# d(profit*)/d{g1_TF24, beta2, leaf_specific_conductance_max} at the optimised
# leaf operating point. Extends test-tf24-leaf-gradient.R (vcmax_25) to the
# hydraulic traits:
#   - g1_TF24, beta2 enter only the hydraulic cost (no transport / assimilation
#     change), so the envelope theorem reduces each to minus the explicit cost
#     derivative (Leaf::dprofit_dg1_TF24 / dprofit_dbeta2);
#   - the supply conductance k_max (= leaf_specific_conductance_max_, which the
#     TF24 trait K_s scales linearly) is a TRANSPORT trait: it moves psi_stem and
#     ci, handled by the transport + IFT pattern (Leaf::dprofit_dkmax).
# Compiled into plant.so, exposed on the Leaf R class -- runs in plain R.

# Standalone leaf at a healthy INTERIOR optimum (g1_TF24 = 5: positive profit,
# dprofit/dcollar ~ 0 where the envelope theorem applies). The same harness as
# test-tf24-leaf-gradient.R; here `run` returns the solved leaf so callers read
# either profit_ or a gradient at the optimised collar.
mk <- function(g1 = 5, beta2 = 1, vc = 100) {
  rc <- 2.65; rb <- 1.29
  Leaf(vcmax_25 = vc, jmax_25 = 167, c = 2.04, b = 3, psi_crit = 5,
       root_c = rc, root_b = rb, root_psi_crit = rb * (log(1 / 0.05))^(1 / rc),
       beta2 = beta2, hk_s = 75, a = 0.3, curv_fact_elec_trans = 0.7,
       curv_fact_colim = 0.99, GSS_tol_abs = 1e-9,
       vulnerability_curve_ncontrol = 100, ci_abs_tol = 1e-6, ci_niter = 1000,
       g1_TF24 = g1, beta_R_H = 3.4e3, beta_R_V = 9.4e4)
}
theta <- 0.000157 * 20; h <- 5
run <- function(l, kmax = theta / h) {
  l$set_physiology(area_leaf = 0.05, mass_root_prop = 1, rho = 608, a_bio = 0.0245,
                   PPFD = 1500, psi_soil = 0.1, soil_depth = 1,
                   leaf_specific_conductance_max = kmax, atm_vpd = 1, ca = 40,
                   sapwood_volume_per_leaf_area = theta * h, leaf_temp = 25,
                   atm_o2_kpa = 21, atm_kpa = 101.3)
  l$find_root_collar_psi()
  l
}

test_that("dprofit_dg1_TF24 matches a re-optimising finite difference", {
  g0 <- 5
  l0 <- run(mk(g1 = g0))
  opt <- -l0$root_collar_psi_
  expect_lt(abs(l0$dprofit_droot_collar_psi(opt)), 1e-3)  # interior optimum
  ad <- l0$dprofit_dg1_TF24(opt)
  profit_at <- function(g) run(mk(g1 = g))$profit_
  hfd <- 1e-5 * g0
  fd <- (profit_at(g0 + hfd) - profit_at(g0 - hfd)) / (2 * hfd)
  expect_equal(ad, fd, tolerance = 1e-5)
})

test_that("dprofit_dbeta2 matches a re-optimising finite difference", {
  b0 <- 1
  l0 <- run(mk(beta2 = b0))
  opt <- -l0$root_collar_psi_
  expect_lt(abs(l0$dprofit_droot_collar_psi(opt)), 1e-3)
  ad <- l0$dprofit_dbeta2(opt)
  profit_at <- function(bb) run(mk(beta2 = bb))$profit_
  hfd <- 1e-5
  fd <- (profit_at(b0 + hfd) - profit_at(b0 - hfd)) / (2 * hfd)
  expect_equal(ad, fd, tolerance = 1e-5)
})

test_that("dprofit_dkmax (transport trait) matches a re-optimising FD", {
  k0 <- theta / h
  l0 <- run(mk(), kmax = k0)
  opt <- -l0$root_collar_psi_
  expect_lt(abs(l0$dprofit_droot_collar_psi(opt)), 1e-3)
  ad <- l0$dprofit_dkmax(opt)
  # re-optimising FD: perturb the supply conductance, re-solve, read profit*.
  profit_at <- function(k) run(mk(), kmax = k)$profit_
  hfd <- 1e-6 * k0
  fd <- (profit_at(k0 + hfd) - profit_at(k0 - hfd)) / (2 * hfd)
  expect_equal(ad, fd, tolerance = 1e-4)
})
