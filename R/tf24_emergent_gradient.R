##' Reverse-mode trait gradient of an SCM's emergent \code{offspring_production}
##' for the TF24 strategy (#472 scope B, Phase F1-full).
##'
##' The TF24 counterpart of \code{\link{offspring_production_gradient}} (FF16). Given a
##' TF24 \code{SCM} run with \code{control(shading_model = "crown-centre",
##' save_RK45_cache = TRUE)}, this returns \eqn{d(\mathrm{offspring\_production}) /
##' d(\theta_k)} for a set of TF24 traits in ONE reverse-mode sweep per cohort -- the
##' calibration-objective gradient (many traits in, one scalar out).
##'
##' TF24 net production comes from a hydraulic leaf optimisation, which has no adjoint
##' tape. So a first (double) pass harvests the trait-independent leaf operating point
##' at every RK stage (the optimised profit, the leaf sensitivities, the
##' \eqn{d(\mathrm{profit})/d(\mathrm{height})} Jacobian); a second pass replays each
##' cohort's survival-weighted offspring over that harvest as a leaf-opt-free, tapeable
##' expression and takes one backward sweep. The full mass cascade + leaf path, the
##' seedling size \code{height_0} (implicit function theorem) and the recruitment
##' filter (establishment, via the seedling net production) are all differentiated; the
##' resident light is held frozen (the rare-mutant / invasion-fitness gradient).
##'
##' The cached resident SCM MUST have been run with \code{shading_model =
##' "crown-centre"} (the replay re-solves the crown-centre leaf optimisation); a
##' deep-crown resident would not be reproduced faithfully.
##'
##' @title Reverse-mode gradient of emergent offspring_production (TF24)
##' @param scm A TF24 \code{SCM} run with \code{control(shading_model = "crown-centre",
##'   save_RK45_cache = TRUE)}. The cached schedule + per-RK-stage resident environment
##'   are read from its patch; the SCM is not re-run.
##' @param traits Character vector of TF24 trait names to differentiate. \code{NULL}
##'   (default) uses all 27 net-production traits (10 leaf + 17 mass-cascade). The
##'   reverse sweep covers all 27 regardless; \code{traits} only selects the output.
##' @param species Integer index of the species (cohort family) to differentiate, in
##'   a multi-species stand. Default \code{1}. \code{offspring_production} is a
##'   per-species emergent output; this returns \eqn{d(\mathrm{offspring\_production}_s)
##'   / d(\theta_k)} for the traits of species \code{s} against the shared frozen
##'   canopy of ALL species (the rare-mutant / invasion gradient).
##' @param birth_rate The (constant) birth-rate driver used in the run. By default it
##'   is recovered as \code{offspring_production / net_reproduction_ratio}.
##' @return A named numeric vector of trait derivatives, with attribute
##'   \code{"offspring_production"} (the value reconstructed by the replay, which should
##'   match \code{scm$offspring_production[[species]]}).
##' @seealso \code{\link{offspring_production_gradient}} (FF16).
##' @export
tf24_offspring_production_gradient <- function(scm, traits = NULL, species = 1L,
                                               birth_rate = NULL) {
  if (is.null(traits)) traits <- tf24_default_traits()
  h <- tf24_harvest(scm, species, birth_rate)
  tf24_offspring_production_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv,
                                          h$ppsab, h$tw, traits)
}

# The 27 net-production TF24 trait names (10 leaf + 17 mass-cascade) the emergent
# gradients differentiate by default.
tf24_default_traits <- function() {
  c("vcmax_25","g1_TF24","beta2","K_s","b","c","jmax_25","a","curv_elec",
    "curv_colim","lma","rho","a_b1","r_l","r_b","r_s","r_r","k_l","k_b",
    "k_s","k_r","a_bio","a_y","a_l1","a_l2","theta","a_r1")
}

# Harvest a run-with-cache TF24 SCM into the frozen pieces the two-pass replay
# consumes (the ResidentHarvest seam, TF24 mirror of ff16_harvest). Requires the
# resident to have been run with shading_model = "crown-centre" (the replay re-solves
# the crown-centre leaf optimisation).
tf24_harvest <- function(scm, species = 1L, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "TF24")) {
    stop("TF24 stand gradients are implemented for the TF24 strategy only")
  }
  if (species < 1L || species > length(scm$patch$species)) {
    stop("species index out of range: stand has ", length(scm$patch$species),
         " species")
  }
  # Cache the patch once: `scm$patch` rebuilds the whole RcppR6 patch object on each
  # access, so the per-stage pr_survival loop below is O(stand size) per call otherwise.
  patch <- scm$patch
  sh <- patch$step_history
  eh <- patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  sp    <- patch$species[[species]]
  nt    <- sp$node_times
  pdens <- sp$patch_densities
  ppsab <- sp$pr_patch_survival_at_birth
  pp    <- unlist(scm$parameters$strategies[[species]]$pars)

  if (is.null(birth_rate)) {
    birth_rate <- scm$offspring_production[[species]] / scm$net_reproduction_ratios[[species]]
  }

  # Cohort birth steps (introductions land on step times).
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
    ppsurv[k, s] <- patch$pr_survival(sh[k] + ah[s] * hN[k])
  }

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, ppsurv = ppsurv,
       ppsab = ppsab, tw = tw, pdens = pdens, nt = nt, birth_rate = birth_rate)
}
