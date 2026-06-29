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
##' @param species Integer index of the species (cohort family) to differentiate, in
##'   a multi-species stand. Default \code{1}. \code{offspring_production} is a
##'   per-species emergent output; this returns \eqn{d(\mathrm{offspring\_production}_s)
##'   / d(\theta_k)} for the traits of species \code{s}. The resident light is the
##'   shared frozen canopy of ALL species (so this is the rare-mutant / invasion
##'   gradient of species \code{s} against the fixed N-species canopy); a
##'   cross-species Jacobian would require a resident-reshaping (active-knot) treatment.
##' @param birth_rate The (constant) birth-rate driver used in the run. By default it
##'   is recovered as \code{offspring_production / net_reproduction_ratio} (exact for
##'   a constant birth rate); pass it explicitly for a time-varying driver.
##' @return A named numeric vector of trait derivatives, with attribute
##'   \code{"offspring_production"} (the value reconstructed by the replay, which
##'   should match \code{scm$offspring_production[[species]]}).
##' @export
offspring_production_gradient <- function(scm, traits = NULL, species = 1L,
                                          birth_rate = NULL) {
  if (is.null(traits)) traits <- ff16_default_traits()
  h <- ff16_harvest(scm, species, birth_rate)
  ff16_offspring_production_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv,
                                          h$ppsab, h$tw, traits)
}

# The 28 production-relevant FF16 trait (parameter) names the emergent gradients
# differentiate by default.
ff16_default_traits <- function() {
  c("lma","rho","theta","a_b1","a_r1","eta_c","a_p1","a_p2","r_l","r_s",
    "r_b","r_r","k_l","k_b","k_s","k_r","a_bio","a_y","a_l1","a_l2",
    "a_f1","a_f2","hmat","omega","a_f3","d_I","a_dG1","a_dG2")
}

# Harvest a run-with-cache FF16 SCM into the frozen pieces the two-pass replay
# consumes (the ResidentHarvest seam): the step schedule, the per-RK-stage resident
# environment, and species `species`'s cohort family (birth steps, node-spacing
# trapezoid weights, per-stage patch survival, survival-at-birth). Shared by every
# FF16 stand-gradient entry point so a new SCM variant only has to emit these.
ff16_harvest <- function(scm, species = 1L, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "FF16")) {
    stop("FF16 stand gradients are implemented for the FF16 strategy only")
  }
  if (species < 1L || species > length(scm$patch$species)) {
    stop("species index out of range: stand has ", length(scm$patch$species),
         " species")
  }
  # Cache the patch once: `scm$patch` rebuilds the whole RcppR6 patch object (every
  # node + species) on each access, so repeatedly indexing it -- especially the
  # per-stage pr_survival loop below -- is O(stand size) per call (~1600x slower).
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
    # Constant birth rate: offspring_production = birth_rate * net_reproduction_ratio.
    birth_rate <- scm$offspring_production[[species]] / scm$net_reproduction_ratios[[species]]
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
    ppsurv[k, s] <- patch$pr_survival(sh[k] + ah[s] * hN[k])
  }

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, ppsurv = ppsurv,
       ppsab = ppsab, tw = tw, pdens = pdens, nt = nt, birth_rate = birth_rate)
}

##' Reverse-mode trait gradient of an SCM's emergent stand metrics (#472 scope B,
##' the calibration-facing generic engine, FF16).
##'
##' The generic counterpart of \code{\link{offspring_production_gradient}}: given a
##' resident \code{SCM} run with \code{control(save_RK45_cache = TRUE)}, it returns a
##' \strong{metrics x traits Jacobian} \eqn{d(\mathrm{metric}_m)/d(\theta_k)} for a
##' set of emergent stand metrics, computed from ONE resident baseline. Every metric
##' is a weighted reduction over the replayed cohorts, \eqn{\mathrm{metric} = \sum_i
##' w_i\, f(\mathrm{state}_i)}; \code{offspring_production} is just one such entry --
##' none is privileged. The engine records one forward replay onto a single adjoint
##' tape and takes one cheap reverse sweep \emph{per metric}, so M metrics cost a
##' replay plus M sweeps, not M replays. This is the calibration core: \code{plant}
##' returns the Jacobian; a downstream package composes likelihoods (data never
##' enters here), which is what lets many likelihood terms share one stand baseline.
##'
##' @title Reverse-mode Jacobian of emergent stand metrics (FF16)
##' @param scm An \code{SCM} run with \code{save_RK45_cache = TRUE} (FF16 strategy).
##' @param metrics Character vector of stand-metric names, any of
##'   \code{"offspring_production"} (the seed-rain integral), \code{"LAI"} (leaf-area
##'   index = the SCM's \code{compute_competition(0)}), \code{"biomass"} (the
##'   size-distribution integral of live + heartwood mass) and \code{"size_moment"}
##'   (the first moment of the size distribution, \eqn{\int n(h)\,h\,dh}). All are
##'   symmetric reductions over the replayed cohorts; none is privileged.
##' @param traits Character vector of FF16 trait names. \code{NULL} (default) uses
##'   all 28 production-relevant parameters.
##' @param species Integer index of the species (cohort family); see
##'   \code{\link{offspring_production_gradient}}. Default \code{1}.
##' @param birth_rate The (constant) birth-rate driver; recovered from the run by
##'   default.
##' @return A list with \code{$jacobian} (a metrics x traits matrix) and
##'   \code{$values} (the reconstructed metric values, which should match the SCM's
##'   emergent outputs).
##' @details Works for both FF16 and TF24 residents (dispatched on the strategy). A
##'   TF24 resident must have been run with \code{shading_model = "crown-centre"}.
##'   FF16 supports all metrics; TF24 currently supports \code{"offspring_production"}
##'   (its census metrics need a leaf-optimisation cross-sensitivity that is a
##'   follow-up). \code{\link{stand_state_jacobian}} works for both.
##' @seealso \code{\link{offspring_production_gradient}}.
##' @export
stand_gradient <- function(scm, metrics = "offspring_production", traits = NULL,
                           species = 1L, birth_rate = NULL) {
  strat <- extract_RcppR6_template_types(scm$parameters, "Parameters")[[1]]
  if (identical(strat, "FF16")) {
    if (is.null(traits)) traits <- ff16_default_traits()
    h <- ff16_harvest(scm, species, birth_rate)
    ff16_stand_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw,
                             traits, metrics, h$birth_rate)
  } else if (identical(strat, "TF24")) {
    if (is.null(traits)) traits <- tf24_default_traits()
    h <- tf24_harvest(scm, species, birth_rate)
    tf24_stand_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw,
                             traits, metrics, h$birth_rate)
  } else {
    stop("stand_gradient is implemented for the FF16 and TF24 strategies only")
  }
}

##' Per-cohort state x trait Jacobian of a resident SCM (#472 scope B, the
##' calibration-facing engine's escape hatch, FF16).
##'
##' The escape hatch for emergent metrics that are NOT a simple weighted reduction
##' (quantiles, ratios, bespoke statistics): rather than register a \code{(w, f)}
##' metric, this exposes \eqn{d(\mathrm{state}_{i,c})/d(\theta_k)} -- the derivative
##' of each replayed cohort's final demographic state component
##' \code{c} \eqn{\in} \{height, mortality, fecundity, area_heartwood,
##' mass_heartwood, offspring\} w.r.t. each trait. ANY smooth downstream metric over
##' the cohort states then has its gradient by the chain rule, with \code{plant}
##' never needing to know the metric -- the same boundary as "likelihoods live
##' downstream". Each cohort's final state is independent, so this tapes one cohort
##' at a time (one reverse sweep per state component).
##'
##' @title Per-cohort state x trait Jacobian (FF16)
##' @param scm An \code{SCM} run with \code{save_RK45_cache = TRUE} (FF16 strategy).
##' @param traits Character vector of FF16 trait names; \code{NULL} uses all 28.
##' @param species Integer species index (see \code{\link{stand_gradient}}).
##' @param birth_rate The (constant) birth-rate driver; recovered by default.
##' @return A list with \code{$states} (a cohort x component matrix of final-state
##'   values) and \code{$jacobian} (a cohort x component x trait array).
##' @details Works for both FF16 and TF24 residents (dispatched on the strategy). A
##'   TF24 resident must have been run with \code{shading_model = "crown-centre"}.
##' @seealso \code{\link{stand_gradient}}.
##' @export
stand_state_jacobian <- function(scm, traits = NULL, species = 1L,
                                 birth_rate = NULL) {
  strat <- extract_RcppR6_template_types(scm$parameters, "Parameters")[[1]]
  if (identical(strat, "FF16")) {
    if (is.null(traits)) traits <- ff16_default_traits()
    h <- ff16_harvest(scm, species, birth_rate)
    ff16_state_jacobian_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw,
                             traits)
  } else if (identical(strat, "TF24")) {
    if (is.null(traits)) traits <- tf24_default_traits()
    h <- tf24_harvest(scm, species, birth_rate)
    tf24_state_jacobian_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw,
                             traits, h$birth_rate)
  } else {
    stop("stand_state_jacobian is implemented for the FF16 and TF24 strategies only")
  }
}
