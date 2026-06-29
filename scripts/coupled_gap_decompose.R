# Step 1 (#472 scope B): decompose the coupled-AD vs full-SCM-FD gap.
#
# The coupled replay freezes THREE channels a full re-run lets move:
#   (a) adaptive ODE step sizes, (b) env-spline knot x-positions, (c) birth-env
#       establishment.
# (a)+(b) are a GRID response: they vanish as the ODE tolerance -> 0 (both the
# resident schedule and the perturbed full runs approach the continuous limit).
# (c) is a real physical channel the replay omits; it does NOT shrink with
# tolerance. So sweeping tolerance separates "grid response" (tol-dependent) from
# "birth-env + boundary" (tol-independent residual).
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

traits <- c("lma", "a_p1", "a_l1")

gap_at_tol <- function(tol_rel, tol_abs) {
  cc <- function(...) control(ode_tol_rel = tol_rel, ode_tol_abs = tol_abs,
                              ode_step_size_min = 1e-9, ...)
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  pr <- run_scm(p, Environment("FF16"), cc(), refine_schedule = TRUE)$parameters
  scm <- run_scm(pr, Environment("FF16"), cc(save_RK45_cache = TRUE),
                 refine_schedule = FALSE)
  h <- plant:::ff16_harvest(scm, 1L, NULL); patch <- scm$patch
  g <- plant:::ff16_coupled_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv,
         h$ppsab, h$tw, traits, "LAI", h$birth_rate,
         patch$stand_newnode_height_stage_history,
         patch$stand_newnode_competition_stage_history, h$patch_area)

  lai_full <- function(trait, delta) {
    q <- pr; q$strategies[[1]]$pars[[trait]] <- q$strategies[[1]]$pars[[trait]] + delta
    s2 <- run_scm(q, Environment("FF16"),
                  control(ode_tol_rel = tol_rel, ode_tol_abs = tol_abs,
                          save_RK45_cache = TRUE), refine_schedule = FALSE)
    e <- s2$patch$environment_history[[length(s2$patch$environment_history)]][[6]]
    -log(e$get_environment_at_height(0))
  }
  out <- data.frame(trait = traits, ad = NA_real_, fd = NA_real_, reldiff = NA_real_)
  for (i in seq_along(traits)) {
    tr <- traits[i]
    d <- 1e-3 * abs(pr$strategies[[1]]$pars[[tr]])
    fd <- (lai_full(tr, +d) - lai_full(tr, -d)) / (2 * d)
    ad <- g$jacobian["LAI", tr]
    out$ad[i] <- ad; out$fd[i] <- fd
    out$reldiff[i] <- abs(ad - fd) / max(abs(fd), 1e-30)
  }
  out
}

cat("=== default tolerance (grid response present) ===\n")
g_def <- gap_at_tol(control()$ode_tol_rel, control()$ode_tol_abs)
print(g_def, row.names = FALSE)

cat("\n=== tight tolerance 1e-7 (grid response suppressed) ===\n")
g_tight <- gap_at_tol(1e-7, 1e-7)
print(g_tight, row.names = FALSE)

cat("\n=== interpretation: residual gap that PERSISTS at tight tol = birth-env+boundary ===\n")
for (i in seq_along(traits)) {
  cat(sprintf("  %-5s  gap(default)=%.2e  gap(tight)=%.2e  grid-response=%.0f%%  birth-env-resid=%.0f%%\n",
    traits[i], g_def$reldiff[i], g_tight$reldiff[i],
    100 * (g_def$reldiff[i] - g_tight$reldiff[i]) / max(g_def$reldiff[i], 1e-30),
    100 * g_tight$reldiff[i] / max(g_def$reldiff[i], 1e-30)))
}
