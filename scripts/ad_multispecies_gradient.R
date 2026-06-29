# Multi-species emergent trait gradient (#472 scope B, Phase C / F1-full).
#
# offspring_production is a PER-SPECIES emergent output. The committed reverse-mode
# APIs -- offspring_production_gradient() (FF16) and tf24_offspring_production_gradient()
# (TF24) -- take a `species` index and return
#
#     d(offspring_production[s]) / d(theta_k)
#
# for the traits theta of species s. The key structural fact that makes the
# single-species machinery generalise by a plain loop: the resident light
# (Patch$environment_history) is the SHARED frozen canopy of ALL species, so every
# species' cohorts replay against the SAME per-RK-stage environment. Holding that
# joint canopy frozen and perturbing only species s's traits is exactly the
# rare-mutant / invasion-fitness gradient of species s against the fixed N-species
# stand -- the natural per-species calibration / selection-gradient target.
#
# (A cross-species Jacobian d(offspring_production[s])/d(theta_{s'}), s' != s, is
# ZERO under the frozen-resident reading -- perturbing another species' traits does
# not move species s's frozen environment. A nonzero cross term only appears once the
# resident canopy is allowed to RESHAPE with the trait, which is the distinct
# active-knot self-shading quantity, e.g. scripts/ad_self_shading_timeint.R.)
#
# Requires `plant` INSTALLED from this branch (run `R CMD INSTALL .` first), plus
# `odelia` (the XAD adjoint tape, resolved at load) and `BH`.
#   Rscript scripts/ad_multispecies_gradient.R

suppressMessages(library(plant))

## ---- A two-species FF16 resident SCM (distinct lma), one shared canopy --------
p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(c(0.0825, 0.2178), "lma"),
                    hyperpar = FF16_hyperpar, birth_rate = list(20, 20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

n_sp <- length(scm$patch$species)
cat(sprintf("Stand: %d species, offspring_production = (%s)\n", n_sp,
            paste(sprintf("%.4f", scm$offspring_production), collapse = ", ")))

## ---- A two-pass FD over species s's frozen schedule (for validation) ----------
# Re-implements the gathering inside offspring_production_gradient() so we can
# finite-difference the reconstructed value, perturbing one trait of species s.
setup <- function(scm, s) {
  sh <- scm$patch$step_history; eh <- scm$patch$environment_history
  sp <- scm$patch$species[[s]]; nt <- sp$node_times
  pp <- unlist(scm$parameters$strategies[[s]]$pars)
  br <- scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]]
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
  N <- length(eh); tcoef <- numeric(length(nt)); x <- nt; n <- length(x)
  tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[n] <- 0.5 * (x[n] - x[n - 1])
  if (n > 2) tcoef[2:(n - 1)] <- 0.5 * (x[3:n] - x[1:(n - 2)])
  tw <- tcoef * sp$patch_densities * pp[["S_D"]] * br
  ah <- c(0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh); ppsurv <- matrix(0, N, 6)
  for (k in seq_len(N)) for (j in 1:6) ppsurv[k, j] <- scm$patch$pr_survival(sh[k] + ah[j] * hN[k])
  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, ppsurv = ppsurv,
       ppsab = sp$pr_patch_survival_at_birth, tw = tw)
}
fd_trait <- function(d, trait, ad, rel_h = c(1e-5, 1e-6)) {
  J_at <- function(v) {
    q <- d$pp; q[[trait]] <- v
    attr(plant:::ff16_offspring_production_gradient_impl(
           q, d$eh, d$sh, d$birth_step, d$ppsurv, d$ppsab, d$tw, trait),
         "offspring_production")
  }
  fds <- vapply(rel_h, function(rh) {
    h <- rh * abs(d$pp[[trait]]); (J_at(d$pp[[trait]] + h) - J_at(d$pp[[trait]] - h)) / (2 * h)
  }, numeric(1))
  fds[which.min(abs(fds - ad))]   # FD truncation/roundoff sweet spot
}

## ---- Per-species gradient: AD (one reverse sweep) vs two-pass FD ---------------
for (s in seq_len(n_sp)) {
  g <- offspring_production_gradient(scm, traits = c("a_p1", "lma"), species = s)
  recon <- attr(g, "offspring_production")
  cat(sprintf("\nSpecies %d (lma = %.4f):\n", s,
              scm$parameters$strategies[[s]]$pars$lma))
  cat(sprintf("  reconstruction: %.6f vs SCM %.6f  (rel.err %.1e)\n",
              recon, scm$offspring_production[[s]],
              abs(recon - scm$offspring_production[[s]]) / scm$offspring_production[[s]]))
  d <- setup(scm, s)
  for (tr in c("a_p1", "lma")) {
    fd <- fd_trait(d, tr, g[[tr]])
    cat(sprintf("  d/d(%-4s): AD %12.5g   FD %12.5g   rel.err %.1e\n",
                tr, g[[tr]], fd, abs(g[[tr]] - fd) / max(1, abs(fd))))
  }
}
