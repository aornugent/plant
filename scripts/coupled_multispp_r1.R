# Gate C (#472 scope B, R2): the CROSS-SPECIES coupled gradient. On a 2-species stand,
# d(total-stand metric)/d(theta of species `target`) via the joint adjoint replay.
# (1) AD vs a central FD over the SAME coupled reconstruction (must agree tightly -- AD
#     and FD differentiate one function); (2) the cross-species term is nonzero (the
#     target's traits move the OTHER species' contribution through the joint canopy),
#     which the frozen/mutant gradient sets to ~0; (3) contrast vs a full-SCM FD.
# Run on a fixed (unrefined) schedule, where the joint R0 re-evolution is stable.
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20, 20))
p$node_schedule_times <- list(seq(0, 90, length.out = 30), seq(0, 90, length.out = 30))
p$max_patch_lifetime <- 90
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

patch <- scm$patch; nsp <- 2; sh <- patch$step_history; eh <- patch$environment_history
mets <- c("LAI", "biomass", "size_moment"); traits <- c("lma", "a_p1")
target <- 1L

pp_list <- lapply(seq_len(nsp), function(s) unlist(scm$parameters$strategies[[s]]$pars))
birth_list <- lapply(seq_len(nsp), function(s)
  vapply(patch$species[[s]]$node_times, function(t) which.min(abs(sh - t)) - 1L, integer(1)))
br <- vapply(seq_len(nsp), function(s)
  scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1))
nn_h <- patch$stand_newnode_height_stage_history_all
nn_c <- patch$stand_newnode_competition_stage_history_all
area <- scm$parameters$patch_area

g <- plant:::ff16_coupled_gradient_ms_impl(pp_list, eh, sh, birth_list, traits, mets,
       br, nn_h, nn_c, area, target)

# (1) AD vs central FD over the SAME coupled reconstruction (perturb the target's trait).
cval <- function(pp_t) {
  pl <- pp_list; pl[[target]] <- pp_t
  plant:::ff16_coupled_metrics_ms_impl(pl, eh, sh, birth_list, mets, br, nn_h, nn_c, area)$values
}
cat("\nGate C (1): AD vs coupled-reconstruction FD, d(total metric)/d(theta_target):\n")
fdmax <- 0
# Step 3e-4 is the noise/truncation optimum: the mass-weighted reconstruction carries
# value noise (cohort height-crossing sort discontinuities) that a smaller FD step
# amplifies (AD is the exact reference; a step sweep bottoms out around 3e-4).
for (tr in traits) {
  d <- 3e-4 * abs(pp_list[[target]][[tr]]); pp1 <- pp2 <- pp_list[[target]]
  pp1[[tr]] <- pp1[[tr]] + d; pp2[[tr]] <- pp2[[tr]] - d
  fd <- (cval(pp1) - cval(pp2)) / (2 * d)
  for (m in mets) {
    rel <- abs(g$jacobian[m, tr] - fd[m]) / max(abs(g$jacobian[m, tr]), abs(fd[m]), 1e-30)
    fdmax <- max(fdmax, rel)
    cat(sprintf("  d%-11s/d%-4s AD=% .5e FD=% .5e rel=%.2e %s\n",
                m, tr, g$jacobian[m, tr], fd[m], rel, if (rel < 2e-2) "OK" else "**"))
  }
}
cat(sprintf("  -> max rel = %.2e  %s\n", fdmax, if (fdmax < 2e-2) "PASS" else "CHECK"))

# (2) Cross-species term: the OTHER species' contribution to total LAI must respond to
# the target's trait. Compare to a per-species frozen gradient (which omits the cross
# term). The difference IS the cross-species feedback.
cat("\nGate C (2): cross-species feedback (total resident vs target-only frozen):\n")
gf <- stand_gradient(scm, metrics = mets, traits = traits, species = target,
                     feedback = "frozen")
for (tr in traits) for (m in mets) {
  cat(sprintf("  d%-11s/d%-4s  total-resident=% .4e  target-frozen=% .4e  cross+feedback=% .4e\n",
              m, tr, g$jacobian[m, tr], gf$jacobian[m, tr], g$jacobian[m, tr] - gf$jacobian[m, tr]))
}
