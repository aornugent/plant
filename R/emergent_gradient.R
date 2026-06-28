##' Reverse-mode trait gradient of an SCM's emergent \code{offspring_production}
##' (#472 scope B, FF16 only).
##'
##' Given an \code{SCM} that has been run with \code{control(save_RK45_cache =
##' TRUE)}, this returns \eqn{d(\mathrm{offspring\_production}) / d(\theta_k)} for a
##' set of FF16 traits \eqn{\theta_k} in a SINGLE reverse-mode sweep -- at the cost
##' of one extra model evaluation, independent of the number of traits, whereas a
##' finite-difference Jacobian needs a fresh whole-stand replay per trait. This is
##' the calibration-objective gradient: many traits in, one scalar out.
##'
##' It is a two-pass method. Pass 1 is the resident SCM run you pass in (it owns the
##' frozen schedule and the per-RK-stage resident light, harvested into
##' \code{patch$step_history} / \code{environment_history}). Pass 2 replays each
##' cohort's demography under the XAD adjoint tape over that frozen schedule
##' (deep-crown assimilation), accumulating the survival-weighted offspring, then
##' takes one backward sweep. The resident light is held frozen (the rare-mutant /
##' invasion-fitness gradient); the recruitment-filter (establishment) initial
##' condition is held frozen too (a separable partial).
##'
##' @title Reverse-mode gradient of emergent offspring_production (FF16)
##' @param scm An \code{SCM} object that has been run with \code{save_RK45_cache =
##'   TRUE} (FF16 strategy). The cached schedule + resident light are read from its
##'   patch; the SCM is not re-run.
##' @param traits Character vector of FF16 trait (parameter) names to differentiate.
##'   \code{NULL} (default) uses all 28 production-relevant parameters.
##' @param birth_rate The (constant) birth-rate driver used in the run. By default it
##'   is recovered as \code{offspring_production / net_reproduction_ratio} (exact for
##'   a constant birth rate); pass it explicitly for a time-varying driver.
##' @return A named numeric vector of trait derivatives, with attribute
##'   \code{"offspring_production"} (the value reconstructed by the replay, which
##'   should match \code{scm$offspring_production}).
##' @export
offspring_production_gradient <- function(scm, traits = NULL, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "FF16")) {
    stop("offspring_production_gradient is implemented for the FF16 strategy only")
  }
  sh <- scm$patch$step_history
  eh <- scm$patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  sp <- scm$patch$species[[1]]
  nt    <- sp$node_times
  pdens <- sp$patch_densities
  ppsab <- sp$pr_patch_survival_at_birth
  pp    <- unlist(scm$parameters$strategies[[1]]$pars)

  if (is.null(birth_rate)) {
    # Constant birth rate: offspring_production = birth_rate * net_reproduction_ratio.
    birth_rate <- scm$offspring_production[[1]] / scm$net_reproduction_ratios[[1]]
  }
  if (is.null(traits)) {
    traits <- c("lma","rho","theta","a_b1","a_r1","eta_c","a_p1","a_p2","r_l","r_s",
                "r_b","r_r","k_l","k_b","k_s","k_r","a_bio","a_y","a_l1","a_l2",
                "a_f1","a_f2","hmat","omega","a_f3","d_I","a_dG1","a_dG2")
  }

  # Cohort birth steps (introductions land exactly on step times).
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))
  N <- length(eh)
  # Node-spacing trapezoid weights so offspring_production == sum_i tw_i * offspring_i.
  tcoef <- numeric(length(nt)); x <- nt; n <- length(x)
  tcoef[1] <- 0.5 * (x[2] - x[1]); tcoef[n] <- 0.5 * (x[n] - x[n - 1])
  if (n > 2) tcoef[2:(n - 1)] <- 0.5 * (x[3:n] - x[1:(n - 2)])
  tw <- tcoef * pdens * pp[["S_D"]] * birth_rate
  # pr_patch_survival at the exact Cash-Karp stage times sh[k] + ah[s]*h.
  ah <- c(0, 0.2, 0.3, 0.6, 1.0, 0.875); hN <- diff(sh)
  ppsurv <- matrix(0, N, 6)
  for (k in seq_len(N)) for (s in 1:6) {
    ppsurv[k, s] <- scm$patch$pr_survival(sh[k] + ah[s] * hN[k])
  }

  ff16_offspring_production_gradient_impl(pp, eh, sh, birth_step, ppsurv, ppsab, tw,
                                          traits)
}
