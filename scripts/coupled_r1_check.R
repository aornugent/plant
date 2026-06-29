# R1 gate for the COUPLED resident replay (#472 scope B). The coupled reverse sweep
# gives the resident TOTAL d(metric)/d(trait) where EVERY trait feeds back (it moves
# heights/densities -> re-shades the active canopy). Validate AD vs a central FD over
# the SAME coupled reconstruction (ff16_coupled_metrics_impl, frozen geometry); AD and
# FD then differentiate one function and must agree tightly. Also contrast with the
# frozen (mutant) gradient -- the resident feedback term, now nonzero for ALL traits
# (vs only a_l1/a_l2 in the frozen-geometry first cut).
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

metrics <- c("LAI", "biomass", "size_moment")
# A spread: allometric (a_l1), a mass/growth trait (lma), an allocation trait (a_p1).
traits  <- c("a_l1", "a_l2", "lma", "a_p1", "hmat")

h <- plant:::ff16_harvest(scm, 1L, NULL)
patch <- scm$patch
nn_h <- patch$stand_newnode_height_stage_history
nn_c <- patch$stand_newnode_competition_stage_history

t0 <- Sys.time()
g <- plant:::ff16_coupled_gradient_impl(
  h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, traits, metrics,
  h$birth_rate, nn_h, nn_c, h$patch_area)
cat(sprintf("coupled reverse sweep: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

coupled_val <- function(pp) {
  plant:::ff16_coupled_metrics_impl(pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab,
    h$tw, metrics, h$birth_rate, nn_h, nn_c, h$patch_area)$values
}
# The coupled reconstruction carries ~1e-8 relative VALUE noise (cohort-height
# crossing sort discontinuities + the env re-evolution drift), so a small FD step is
# noise-dominated; a step sweep shows AD vs FD bottoms out at ~1e-4..1e-6 around
# step 1e-4 (the noise/truncation optimum). AD is the exact derivative; FD is the
# noise-limited reference, so we validate at the optimal step.
fd_col <- function(trait, rel = 1e-4) {
  d <- rel * abs(h$pp[[trait]]); if (d == 0) d <- rel
  pp_p <- h$pp; pp_p[[trait]] <- pp_p[[trait]] + d
  pp_m <- h$pp; pp_m[[trait]] <- pp_m[[trait]] - d
  (coupled_val(pp_p) - coupled_val(pp_m)) / (2 * d)
}

cat("\nR1 AD vs coupled-reconstruction FD (d(metric)/d(trait)):\n")
fdmax <- 0
for (tr in traits) {
  ad <- g$jacobian[, tr]; fd <- fd_col(tr)
  for (m in metrics) {
    rel <- abs(ad[m] - fd[m]) / pmax(abs(ad[m]), abs(fd[m]), 1e-30)
    fdmax <- max(fdmax, rel)
    cat(sprintf("  d%-12s/d%-4s  AD=% .6e  FD=% .6e  rel=%.2e %s\n",
                m, tr, ad[m], fd[m], rel, if (rel < 5e-3) "OK" else "**"))
  }
}

# Contrast with the frozen (mutant) gradient -- the resident feedback term.
gf <- stand_gradient(scm, metrics = metrics, traits = traits, feedback = "frozen")
cat("\nR1 resident TOTAL (coupled) vs frozen (mutant) -- feedback per trait:\n")
for (tr in traits) for (m in metrics) {
  fr <- gf$jacobian[m, tr]; rs <- g$jacobian[m, tr]
  flip <- if (sign(fr) != sign(rs) && fr != 0) "  <-- SIGN FLIP" else ""
  cat(sprintf("  d%-12s/d%-4s  frozen=% .4e  resident=% .4e  feedback=% .4e%s\n",
              m, tr, fr, rs, rs - fr, flip))
}
cat(sprintf("\nR1 AD-vs-FD max rel error = %.2e  %s\n", fdmax,
            if (fdmax < 5e-3) "PASS" else "** CHECK **"))
