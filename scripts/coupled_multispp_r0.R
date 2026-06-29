# Gate B (#472 scope B, R2): multi-species coupled R0. Does the JOINT double-precision
# re-evolution (all species' cohorts stepped together, joint canopy reconstructed each
# RK stage) reproduce (a) the SCM's per-stage env (env_err small) and (b) the TOTAL
# stand metrics (vs the single-species frozen engine summed over species, and LAI vs
# the SCM's compute_competition(0))?
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20, 20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

patch <- scm$patch
nsp <- length(scm$parameters$strategies)
sh  <- patch$step_history
eh  <- patch$environment_history
metrics <- c("LAI", "biomass", "size_moment")

# Per-species harvest.
pp_list <- lapply(seq_len(nsp), function(s)
  unlist(scm$parameters$strategies[[s]]$pars))
birth_list <- lapply(seq_len(nsp), function(s) {
  nt <- patch$species[[s]]$node_times
  vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
})
birth_rate <- vapply(seq_len(nsp), function(s)
  scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1))
nn_h <- patch$stand_newnode_height_stage_history_all
nn_c <- patch$stand_newnode_competition_stage_history_all
area <- scm$parameters$patch_area

ms <- plant:::ff16_coupled_metrics_ms_impl(pp_list, eh, sh, birth_list, metrics,
        birth_rate, nn_h, nn_c, area)

# Cross-check 1: total LAI == SCM compute_competition(0) at the final census.
scm_lai <- -log(eh[[length(eh)]][[6]]$get_environment_at_height(0))

# Cross-check 2: total biomass / size_moment via the single-species FROZEN engine
# (validated vs the SCM) summed over species. The frozen engine reads the joint env, so
# its per-species census is correct; at theta0 the coupled value should match the sum.
frz_tot <- setNames(numeric(length(metrics)), metrics)
for (s in seq_len(nsp)) {
  gf <- stand_gradient(scm, metrics = metrics, species = s, feedback = "frozen")
  frz_tot <- frz_tot + gf$values
}

cat(sprintf("\nMS R0 env drift: env_err = %.3e (z=%.3f)\n", ms$env_err, ms$env_err_z))
cat("\nMS R0 TOTAL metric values (coupled-MS vs frozen-engine-summed):\n")
for (m in metrics) {
  rel <- abs(ms$values[m] - frz_tot[m]) / max(abs(frz_tot[m]), 1e-30)
  cat(sprintf("  %-12s coupled=% .8e  frozen-sum=% .8e  rel=%.2e %s\n",
              m, ms$values[m], frz_tot[m], rel, if (rel < 1e-4) "OK" else "**"))
}
cat(sprintf("  LAI vs SCM competition(0): coupled=% .8e SCM=% .8e rel=%.2e\n",
            ms$values["LAI"], scm_lai, abs(ms$values["LAI"] - scm_lai) / abs(scm_lai)))
