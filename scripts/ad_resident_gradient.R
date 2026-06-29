# Resident TOTAL gradient of emergent stand metrics (#472 scope B, R0-R1).
#
# The shipped stand_gradient() defaults to feedback="frozen": the resident canopy is
# held fixed and each cohort reads the harvested env as a constant (the rare-mutant /
# invasion-fitness gradient). For LAI / biomass / size-moment as RESIDENT ecosystem
# outcomes the canopy co-varies with the trait: an allometric trait (a_l1, a_l2)
# reshapes EVERY resident's leaf area, hence the Beer's-law light every plant reads.
# feedback="resident" turns that channel on via a VALUE-ANCHORED reconstruction over
# the new per-RK-stage stand harvest: each cohort reads the exact frozen-env VALUE
# plus the trait-DERIVATIVE of a trapezium reconstruction (zero value), so the metric
# VALUES are bit-identical to the frozen engine (R0 baseline gate) while the gradient
# gains the resident feedback (R1).
#
# Validations:
#   R0  feedback="resident" $values == feedback="frozen" $values  (anchoring, ~0).
#   R0  per-RK-stage harvest aligns 1:1 with environment_history.
#   R1  AD vs reconstruction-FD on the (un-anchored) recon  -> machinery exact.
#   R1  resident vs frozen gradient -> the committed feedback term (C-27 sign flip).
#
# Run after `make full_compile`:  Rscript scripts/ad_resident_gradient.R
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

metrics <- c("LAI", "biomass", "size_moment")
traits  <- c("a_l1", "a_l2")

## ---- R0: per-RK-stage stand harvest aligns with environment_history -----------
patch <- scm$patch
eh  <- patch$environment_history
shs <- patch$stand_height_stage_history
scs <- patch$stand_competition_stage_history
cat(sprintf("R0 harvest: %d steps; stand-stage history %d steps; per-step stages %s\n",
            length(eh), length(shs),
            paste(unique(vapply(shs, length, integer(1))), collapse = ",")))
stopifnot(length(shs) == length(eh), length(scs) == length(eh))
stopifnot(all(vapply(shs, length, integer(1)) == 6L))
# last-stage stand of step n ~ step-boundary stand (stand_height_history[n]).
sh_step <- patch$stand_height_history
n <- length(eh)
align_err <- max(abs(sort(shs[[n]][[6]]) - sort(sh_step[[n]])))
cat(sprintf("R0 harvest: |last-stage stand - step-boundary stand| (step %d) = %.2e\n",
            n, align_err))

## ---- R0: baseline gate -- resident VALUES == frozen VALUES --------------------
g_froz <- stand_gradient(scm, metrics = metrics, traits = traits, feedback = "frozen")
g_res  <- stand_gradient(scm, metrics = metrics, traits = traits, feedback = "resident")
val_err <- max(abs(g_res$values - g_froz$values))
cat("\nR0 baseline gate (metric VALUES, resident vs frozen):\n")
print(rbind(frozen = g_froz$values, resident = g_res$values))
cat(sprintf("  max |resident - frozen| value = %.2e  %s\n", val_err,
            if (val_err < 1e-10) "OK (value-anchored)" else "** DRIFT **"))

## ---- R1: AD vs reconstruction-FD (un-anchored), machinery exactness -----------
h  <- plant:::ff16_harvest(scm, 1L, NULL)
a_l1_0 <- h$pp[["a_l1"]]; a_l2_0 <- h$pp[["a_l2"]]
call_impl <- function(pp, feedback) {
  plant:::ff16_stand_gradient_impl(pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab,
    h$tw, traits, metrics, h$birth_rate, feedback, h$sh_h, h$sh_c, h$patch_area,
    a_l1_0, a_l2_0)                       # weight basis FROZEN at the resident base
}
g_na <- call_impl(h$pp, "resident_noanchor")   # AD total gradient (genuine recon value)

fd_col <- function(trait, rel = 1e-6) {
  d <- rel * h$pp[[trait]]
  pp_p <- h$pp; pp_p[[trait]] <- pp_p[[trait]] + d
  pp_m <- h$pp; pp_m[[trait]] <- pp_m[[trait]] - d
  (call_impl(pp_p, "resident_noanchor")$values -
   call_impl(pp_m, "resident_noanchor")$values) / (2 * d)
}
cat("\nR1 AD vs reconstruction-FD (un-anchored recon; d(metric)/d(trait)):\n")
for (tr in traits) {
  ad <- g_na$jacobian[, tr]
  fd <- fd_col(tr)
  rel <- abs(ad - fd) / pmax(abs(ad), 1e-30)
  for (m in metrics)
    cat(sprintf("  d%-12s/d%-4s  AD=% .6e  FD=% .6e  rel.err=%.2e %s\n",
                m, tr, ad[m], fd[m], rel[m], if (rel[m] < 1e-4) "OK" else "**"))
}

## ---- R1: the committed resident feedback (resident vs frozen gradient) --------
cat("\nR1 resident TOTAL vs frozen (mutant) gradient -- the resident feedback term:\n")
for (tr in traits) for (m in metrics) {
  fr <- g_froz$jacobian[m, tr]; rs <- g_res$jacobian[m, tr]
  flip <- if (sign(fr) != sign(rs) && fr != 0) "  <-- SIGN FLIP" else ""
  cat(sprintf("  d%-12s/d%-4s  frozen=% .6e  resident=% .6e  feedback=% .6e%s\n",
              m, tr, fr, rs, rs - fr, flip))
}

# z-channel: the shipped "resident" path carries the focal-height->light derivative via
# the frozen env's analytic z-derivative; "resident_noanchor" uses the recon's own
# z-derivative (the FD target). They agree to ~3% (the recon-vs-spline z-derivative gap).
zchan_err <- max(abs(g_res$jacobian[, traits] - g_na$jacobian[, traits]) /
                 pmax(abs(g_na$jacobian[, traits]), 1e-30))
cat(sprintf("\nz-channel (anchored frozen-ld vs noanchor recon-z): max rel = %.2e\n",
            zchan_err))

# R1 gate: the forward resident machinery reproduces the reconstruction-FD tightly.
fd_max <- 0
for (tr in traits) fd_max <- max(fd_max, abs(g_na$jacobian[, tr] - fd_col(tr)) /
                                 pmax(abs(g_na$jacobian[, tr]), 1e-30))
stopifnot(val_err < 1e-10, fd_max < 1e-3)
cat(sprintf("\nR0-R1 resident-gradient validation complete (R1 AD-vs-FD max rel = %.2e).\n",
            fd_max))
