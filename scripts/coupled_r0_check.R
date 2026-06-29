# R0 gate for the COUPLED resident replay (#472 scope B, course-corrected build).
#
# The coupled replay re-evolves the WHOLE stand together over the frozen schedule,
# reconstructing the active canopy light each RK stage from the current active stand
# (heights AND densities respond to theta). R0 asks: at theta0 does this double-
# precision re-evolution reproduce (a) the SCM's per-stage env (env_err ~ 3e-14) and
# (b) the emergent metric VALUES (vs the frozen engine and the SCM itself)?
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

metrics <- c("LAI", "biomass", "size_moment")
h <- plant:::ff16_harvest(scm, 1L, NULL)
patch <- scm$patch
nn_h <- patch$stand_newnode_height_stage_history
nn_c <- patch$stand_newnode_competition_stage_history

cpl <- plant:::ff16_coupled_metrics_impl(
  h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, metrics, h$birth_rate,
  nn_h, nn_c, h$patch_area)

# Frozen engine values (per-cohort independent replay over the frozen env).
frz <- stand_gradient(scm, metrics = metrics, feedback = "frozen")$values

# SCM ground truth: LAI metric == competition(0) at the final census == -log light(0)
# from the final cached env; the engine's other metrics have no direct SCM scalar, so
# we lean on the frozen-engine cross-check (it is validated against the SCM elsewhere).
e_fin <- h$eh[[length(h$eh)]][[6]]
scm_lai <- -log(e_fin$get_environment_at_height(0))

cat(sprintf("\nR0 coupled env drift: env_err = %.3e  %s\n", cpl$env_err,
            if (cpl$env_err < 1e-12) "PASS (re-evolution reproduces SCM stand)"
            else "** coupled drift **"))

cat("\nR0 metric VALUES (coupled vs frozen-engine vs SCM):\n")
for (m in metrics) {
  rel <- abs(cpl$values[m] - frz[m]) / pmax(abs(frz[m]), 1e-30)
  cat(sprintf("  %-12s coupled=% .8e  frozen=% .8e  rel=%.2e %s\n",
              m, cpl$values[m], frz[m], rel, if (rel < 1e-9) "OK" else "**"))
}
cat(sprintf("  LAI vs SCM competition(0): coupled=% .8e  SCM=% .8e  rel=%.2e\n",
            cpl$values["LAI"], scm_lai,
            abs(cpl$values["LAI"] - scm_lai) / abs(scm_lai)))
