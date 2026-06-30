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

  # Per-RK-stage resident stand (species 0) for the RESIDENT feedback path (R0):
  # [step][stage 0..5][cohort] heights + per-node competition effects, aligned 1:1
  # with environment_history. Empty unless the run cached them (older caches lack it).
  sh_h <- patch$stand_height_stage_history
  sh_c <- patch$stand_competition_stage_history
  # Boundary new_node (height + competition effect) per RK stage: the trapezium tail
  # term Species::compute_competition adds beyond `nodes`, needed by the COUPLED
  # resident replay's per-stage canopy reconstruction at ground level.
  nn_h <- patch$stand_newnode_height_stage_history
  nn_c <- patch$stand_newnode_competition_stage_history
  patch_area <- scm$parameters$patch_area

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, ppsurv = ppsurv,
       ppsab = ppsab, tw = tw, pdens = pdens, nt = nt, birth_rate = birth_rate,
       sh_h = sh_h, sh_c = sh_c, nn_h = nn_h, nn_c = nn_c, patch_area = patch_area)
}

# Harvest a run-with-cache FF16 SCM into the ALL-SPECIES arrays the multi-species
# coupled engine consumes (the cross-species resident Jacobian, #472 scope B R2): the
# shared schedule + joint env, per-species parameter vectors / cohort birth steps /
# constant birth rates, and the all-species per-RK-stage boundary harvest
# (stand_newnode_*_stage_history_all, [step][stage][species]). The joint stand light is
# reconstructed in C++ from each species' re-evolved cohorts; only these per-species
# pieces are needed from R.
ff16_harvest_ms <- function(scm) {
  patch <- scm$patch
  nsp <- length(scm$parameters$strategies)
  sh  <- patch$step_history
  eh  <- patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  nn_h <- patch$stand_newnode_height_stage_history_all
  nn_c <- patch$stand_newnode_competition_stage_history_all
  if (length(nn_h) < 1L) {
    stop("feedback = 'resident' on a multi-species stand needs the all-species ",
         "per-RK-stage harvest; re-run the resident SCM on this plant version with ",
         "control(save_RK45_cache = TRUE)")
  }
  list(
    pp_list = lapply(seq_len(nsp), function(s)
      unlist(scm$parameters$strategies[[s]]$pars)),
    eh = eh, sh = sh,
    birth_list = lapply(seq_len(nsp), function(s)
      vapply(patch$species[[s]]$node_times,
             function(t) which.min(abs(sh - t)) - 1L, integer(1))),
    birth_rate = vapply(seq_len(nsp), function(s)
      scm$offspring_production[[s]] / scm$net_reproduction_ratios[[s]], numeric(1)),
    nn_h = nn_h, nn_c = nn_c,
    patch_area = scm$parameters$patch_area, nsp = nsp)
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
##' @param feedback How the resident light responds to the trait. \code{"frozen"}
##'   (default) holds the canopy fixed -- the rare-mutant / invasion-fitness gradient,
##'   correct for \code{offspring_production}. \code{"resident"} is the resident TOTAL
##'   gradient via the COUPLED whole-stand replay: all cohorts are re-evolved together
##'   over the frozen schedule and the canopy light is reconstructed each RK stage from
##'   the active stand (cohort heights AND densities respond to the trait, odelia #32
##'   active-knot spline). So EVERY trait feeds back -- a trait that changes growth or
##'   mortality moves the canopy everyone reads -- and the feedback routinely dominates
##'   and flips the sign of the census metrics (LAI / biomass / size-moment) relative
##'   to the frozen reading. Applies to the census metrics; \code{offspring_production}
##'   stays \code{"frozen"} (the invasion gradient) even under \code{"resident"}. FF16
##'   only so far (#472 scope B). \strong{Recommended use} (scope decision): for the
##'   census metrics (LAI / biomass / size_moment) \code{"resident"} is the correct
##'   stand-level total gradient and \code{"frozen"} the rare-mutant invasion gradient;
##'   for \code{offspring_production} the two coincide (it is always the invasion
##'   gradient). The default is left at \code{"frozen"} so the function is a safe,
##'   backward-compatible no-feedback derivative unless feedback is asked for.
##'   \strong{Multi-species:} on a stand with more than one species \code{"resident"}
##'   returns the CROSS-SPECIES total gradient \eqn{d(\mathrm{total\ stand\ metric})/
##'   d(\theta_{\mathrm{species}})} -- the joint canopy is rebuilt from every species'
##'   re-evolved cohorts and the census metrics are summed over species, so the
##'   differentiated species' traits feed back through the canopy that every species
##'   reads (the cross term \code{"frozen"} sets to zero). The joint re-evolution is
##'   well-conditioned on a \emph{fixed} node schedule; on a strongly clustered
##'   adaptively-refined schedule it can go stiff, and the function then raises a clear
##'   error (gated by a cheap baseline pass) asking for a fixed/uniform schedule rather
##'   than returning a diverged result.
##' @return A list with \code{$jacobian} (a metrics x traits matrix) and
##'   \code{$values} (the reconstructed metric values, which should match the SCM's
##'   emergent outputs).
##' @details Works for FF16, TF24 and TF24f residents (dispatched on the strategy). A
##'   TF24/TF24f resident must have been run with \code{shading_model = "crown-centre"}.
##'   FF16 supports all metrics (frozen + resident). TF24 currently supports
##'   \code{"offspring_production"} (its census metrics need a leaf-optimisation
##'   cross-sensitivity that is a follow-up). TF24f -- the fast-acclimation variant
##'   whose tracked-collar leaf eval is analytic -- supports the CENSUS metrics
##'   (\code{"LAI"} / \code{"biomass"} / \code{"size_moment"}) under \code{feedback =
##'   "frozen"} via the reverse-mode AD census tape (the collar is carried as a taped
##'   state with a curvature-linearised gradient-ascent rate); its resident (coupled)
##'   census gradient and an offspring tape are follow-ups. \code{\link{stand_state_jacobian}}
##'   works for FF16 and TF24.
##' @seealso \code{\link{offspring_production_gradient}}.
##' @export
stand_gradient <- function(scm, metrics = "offspring_production", traits = NULL,
                           species = 1L, birth_rate = NULL,
                           feedback = c("frozen", "resident")) {
  # "resident_noanchor" is an undocumented validation mode (genuine recon value, for
  # FD-checking R1); the public choices are "frozen"/"resident".
  feedback <- if (length(feedback) > 1L) match.arg(feedback) else
    match.arg(feedback, c("frozen", "resident", "resident_noanchor"))
  is_resident <- feedback %in% c("resident", "resident_noanchor")
  strat <- extract_RcppR6_template_types(scm$parameters, "Parameters")[[1]]
  if (identical(strat, "FF16")) {
    if (is.null(traits)) traits <- ff16_default_traits()
    if (feedback == "frozen") {
      # Frozen invasion gradient -> the FULLY native entry: the whole harvest (env +
      # schedule + birth steps + weights + per-stage survival) is built in C++ from the
      # live Patch, so the O(stand) R-side ff16_harvest (pr_survival loop) never runs.
      pp <- unlist(scm$parameters$strategies[[species]]$pars)
      return(ff16_stand_gradient_native(scm, pp, as.integer(species - 1L), traits,
               metrics, if (is.null(birth_rate)) -1 else birth_rate, "frozen",
               list(), list(), scm$parameters$patch_area, -1, -1))
    }
    if (feedback != "resident") {
      # The undocumented resident_noanchor validation mode -> the per-cohort frozen-
      # canopy / leaf-area graft engine (needs the per-RK-stage stand harvest, so this
      # rare path still uses the R harvest).
      h <- ff16_harvest(scm, species, birth_rate)
      if (is_resident && length(h$sh_h) < 1L) {
        stop("feedback = 'resident' needs the per-RK-stage stand harvest; re-run the ",
             "resident SCM on this plant version with control(save_RK45_cache = TRUE)")
      }
      return(ff16_stand_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv,
               h$ppsab, h$tw, traits, metrics, h$birth_rate, feedback,
               if (is.null(h$sh_h)) list() else h$sh_h,
               if (is.null(h$sh_c)) list() else h$sh_c, h$patch_area, -1, -1))
    }
    # feedback = "resident": the COUPLED whole-stand replay (every trait re-shades the
    # canopy). Census metrics (LAI / biomass / size_moment) use the coupled total
    # gradient; offspring_production stays the FROZEN invasion gradient (the canopy a
    # rare mutant invades is the resident's, not co-moving with the mutant's trait).
    # Fully native for single-species: cheap params from $parameters (no scm$patch
    # rebuild, no ff16_harvest); the boundary-node guard lives in the C++ entry.
    pp <- unlist(scm$parameters$strategies[[species]]$pars)
    patch_area <- scm$parameters$patch_area
    br <- if (is.null(birth_rate)) -1 else birth_rate
    census_set <- c("LAI", "biomass", "size_moment")
    bad <- setdiff(metrics, c(census_set, "offspring_production"))
    if (length(bad)) stop("unknown stand metric: ", paste(bad, collapse = ", "))
    census <- metrics[metrics %in% census_set]
    offsp  <- metrics[metrics == "offspring_production"]
    jac <- matrix(0, length(metrics), length(traits),
                  dimnames = list(metrics, traits))
    values <- stats::setNames(numeric(length(metrics)), metrics)
    nsp <- length(scm$parameters$strategies)
    if (length(census)) {
      if (nsp == 1L) {
        # Single-species coupled total gradient (canopy = this species' re-evolved stand).
        # Fully native: env + harvest + boundary-node history from the live Patch.
        gc <- ff16_coupled_gradient_native(scm, pp, as.integer(species - 1L),
                traits, census, br, patch_area)
      } else {
        # Multi-species CROSS-SPECIES coupled Jacobian: the joint canopy is rebuilt from
        # ALL species' re-evolved cohorts and the census metrics are TOTAL-stand sums, so
        # this returns d(total-stand metric)/d(theta of species `species`) -- including
        # the cross term whereby the differentiated species re-shades the canopy every
        # species reads. The joint re-evolution is well-conditioned on a fixed node
        # schedule but can go stiff (one cohort's log-density runs away) on a strongly
        # clustered ADAPTIVELY-REFINED schedule; a cheap double R0 pass gates it and
        # raises a clear error (rather than returning NaN) so the caller can re-run the
        # resident SCM with a fixed/uniform node schedule.
        hm <- ff16_harvest_ms(scm)
        r0 <- ff16_coupled_metrics_ms_impl(hm$pp_list, hm$eh, hm$sh, hm$birth_list,
                census, hm$birth_rate, hm$nn_h, hm$nn_c, hm$patch_area)
        if (!all(is.finite(r0$values)) || r0$env_err > 1e-2) {
          stop("the multi-species coupled resident re-evolution diverged on this node ",
               "schedule (joint env drift = ", signif(r0$env_err, 3), "). The ",
               "cross-species coupled gradient needs a well-conditioned node schedule; ",
               "re-run the resident SCM with a fixed/uniform schedule (e.g. ",
               "p$node_schedule_times <- list(seq(0, T, length.out = n), ...); ",
               "run_scm(p, ..., refine_schedule = FALSE)) rather than an adaptively ",
               "refined one.")
        }
        gc <- ff16_coupled_gradient_ms_impl(hm$pp_list, hm$eh, hm$sh, hm$birth_list,
                traits, census, hm$birth_rate, hm$nn_h, hm$nn_c, hm$patch_area,
                as.integer(species))
      }
      jac[census, ] <- gc$jacobian[census, , drop = FALSE]
      values[census] <- gc$values[census]
    }
    if (length(offsp)) {
      # offspring stays the FROZEN invasion gradient (fully native, like the frozen path).
      go <- ff16_stand_gradient_native(scm, pp, as.integer(species - 1L), traits,
              offsp, br, "frozen", list(), list(), patch_area, -1, -1)
      jac[offsp, ] <- go$jacobian[offsp, , drop = FALSE]
      values[offsp] <- go$values[offsp]
    }
    list(jacobian = jac, values = values)
  } else if (identical(strat, "TF24")) {
    if (is_resident) {
      stop("feedback = 'resident' is implemented for FF16 only so far (R0-R1); ",
           "TF24 resident light is R2")
    }
    if (is.null(traits)) traits <- tf24_default_traits()
    h <- tf24_harvest(scm, species, birth_rate)
    tf24_stand_gradient_impl(h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw,
                             traits, metrics, h$birth_rate)
  } else if (identical(strat, "TF24f")) {
    # TF24f: the CENSUS metrics (LAI / biomass / size_moment) via the reverse-mode
    # AD tape (#472 scope B). The tracked-collar leaf eval is analytic, so -- unlike
    # TF24 -- the census number-density gradient is available (the collar is carried as
    # a taped state with a curvature-linearised rate). feedback = "frozen" gives the
    # rare-mutant / invasion census gradient (tf24f_census_gradient_ad, step 2);
    # feedback = "resident" gives the coupled TOTAL stand gradient where every trait
    # re-shades the canopy (tf24f_resident_census_gradient_ad, step 5). offspring_production
    # (its own tape) is a separate follow-up.
    census_set <- c("LAI", "biomass", "size_moment")
    bad <- setdiff(metrics, c(census_set, "offspring_production"))
    if (length(bad)) {
      stop("stand_gradient for TF24f supports offspring_production and the census ",
           "metrics (", paste(census_set, collapse = ", "), "). Got: ",
           paste(bad, collapse = ", "))
    }
    if (is.null(traits)) traits <- tf24_default_traits()
    census <- metrics[metrics %in% census_set]
    offsp  <- metrics[metrics == "offspring_production"]
    jac <- matrix(0, length(metrics), length(traits),
                  dimnames = list(metrics, traits))
    values <- stats::setNames(numeric(length(metrics)), metrics)
    if (length(census)) {
      # feedback = "frozen" -> the rare-mutant / invasion census gradient (step 2);
      # feedback = "resident" -> the coupled TOTAL stand gradient (step 5, single-species).
      gc <- if (is_resident) {
        if (length(scm$patch$species) > 1L) {
          # Multi-species: the CROSS-SPECIES total gradient d(total-stand metric)/d(theta
          # of species `species`) -- the joint canopy is rebuilt from every species'
          # re-evolved cohorts, so the differentiated species feeds back through the canopy
          # every species reads. Well-conditioned on a fixed node schedule (gated for stiff
          # schedules inside the ms wrapper).
          tf24f_resident_census_gradient_ms_ad(scm, metrics = census, traits = traits,
                                               species = species)
        } else {
          tf24f_resident_census_gradient_ad(scm, metrics = census, traits = traits,
                                            species = species, birth_rate = birth_rate)
        }
      } else {
        tf24f_census_gradient_ad(scm, metrics = census, traits = traits,
                                 species = species, birth_rate = birth_rate)
      }
      jac[census, ] <- gc$jacobian[census, , drop = FALSE]
      values[census] <- gc$values[census]
    }
    if (length(offsp)) {
      # offspring_production is always the FROZEN invasion gradient (the canopy a rare
      # mutant invades is the resident's), even under feedback = "resident".
      go <- tf24f_offspring_production_gradient(scm, traits = traits, species = species,
                                                birth_rate = birth_rate)
      jac["offspring_production", ] <- go$gradient
      values["offspring_production"] <- go$value
    }
    list(jacobian = jac, values = values)
  } else {
    stop("stand_gradient is implemented for the FF16, TF24 and TF24f strategies only")
  }
}

##' Reverse-mode trait gradient of \code{\link{grow_individual_to_size}} (#472 scope
##' B, FF16): a single plant grown in a FIXED environment to target size(s),
##' differentiated w.r.t. traits.
##'
##' The per-plant, fixed-environment counterpart of \code{\link{stand_gradient}}: there
##' is no resident feedback (the environment is given), so the gradient is the exact
##' derivative of the grow-to-size solve. For each target size it returns the derivative
##' of the stopping TIME \eqn{t^*} (the time the plant reaches the target) and of every
##' ODE STATE component at \eqn{t^*}, w.r.t. all FF16 traits, as an optional output of
##' the same call you already make. This is the gradient a growth-rate or
##' time-to-size calibration / optimisation consumes (e.g.
##' \code{\link{optimise_individual_rate_at_size_by_trait}}).
##'
##' It is a two-pass method mirroring the SCM gradients. Pass 1 runs the ordinary
##' \code{\link{grow_individual_to_size}} (its adaptive Cash-Karp step schedule and
##' per-node trajectory are harvested). Pass 2 replays the demographic ODE over that
##' FROZEN schedule with the trait active, reading the fixed environment with the
##' default deep-crown assimilation, to a single partial final step landing on
##' \eqn{t^*}; one reverse sweep per state component gives \eqn{\partial
##' \mathrm{state}/\partial\theta} at fixed \eqn{t^*}. The stopping time itself responds
##' to the trait through the implicit function theorem on the stopping condition
##' \eqn{\mathrm{size}(t^*,\theta) = \mathrm{target}}:
##' \deqn{dt^*/d\theta = -\,(\partial\,\mathrm{size}/\partial\theta\,|_{t^*}) /
##'   \dot{\mathrm{size}}(t^*),}
##' so the TOTAL derivative of each returned component \eqn{y_c} is \eqn{dy_c/d\theta =
##' \partial y_c/\partial\theta|_{t^*} + \dot y_c(t^*)\,dt^*/d\theta} (for the size
##' component itself the two terms cancel, as it is pinned to the target). The seedling
##' size \eqn{h_0} (which solves \eqn{\mathrm{mass}(h_0) = \mathrm{seed\ mass}}) carries
##' its own \eqn{dh_0/d\theta} by the same implicit-function step.
##'
##' @title Reverse-mode gradient of grow_individual_to_size (FF16 / TF24f)
##' @param individual An \code{Individual} object (FF16 or TF24f strategy), as passed to
##'   \code{\link{grow_individual_to_size}}. For TF24f the trajectory carries the tracked
##'   collar (\code{opt_root_psi_state}) as a 6th state and the gradient is dispatched to
##'   the TF24f AD tape (same frozen-schedule replay + stopping-time IFT, with the collar
##'   curvature-linearised); all other behaviour is identical.
##' @param sizes A vector of target sizes to grow the plant to (increasing).
##' @param size_name The size variable the targets refer to (one of the ODE state
##'   names, e.g. \code{"height"}; FF16's monotonic size).
##' @param env An \code{Environment} object (the fixed environment).
##' @param traits Character vector of FF16 trait names; \code{NULL} (default) uses all
##'   28 production-relevant parameters.
##' @param time_max,warn Passed through to \code{\link{grow_individual_to_size}} for the
##'   schedule-discovery pass.
##' @return A list with \code{$time} (the reconstructed \eqn{t^*} per size),
##'   \code{$state} (a sizes x component matrix of the ODE state at \eqn{t^*}),
##'   \code{$d_time} (a sizes x trait matrix \eqn{dt^*/d\theta}) and \code{$d_state} (a
##'   sizes x component x trait array of TOTAL \eqn{d\,\mathrm{state}/d\theta}). The
##'   reconstructed \code{$time}/\code{$state} match \code{grow_individual_to_size} to
##'   the live \code{uniroot} tolerance.
##' @seealso \code{\link{grow_individual_to_size}}, \code{\link{stand_gradient}}.
##' @export
grow_individual_to_size_gradient <- function(individual, sizes, size_name, env,
                                             traits = NULL, time_max = Inf,
                                             warn = TRUE) {
  # TF24f: the fast-acclimation variant carries the tracked collar as a 6th ODE state;
  # its grow gradient reuses the same frozen-schedule replay + stopping-time IFT but with
  # the curvature-linearised collar (build-order step 4). Dispatch to the TF24f AD tape.
  if (grepl("^TF24f", individual$strategy_name)) {
    return(tf24f_grow_individual_to_size_gradient_ad(individual, sizes, size_name, env,
                                                     traits, time_max, warn))
  }
  if (!grepl("^FF16", individual$strategy_name))
    stop("grow_individual_to_size_gradient is implemented for the FF16 and TF24f ",
         "strategies only")
  if (is.unsorted(sizes) || length(sizes) == 0L)
    stop("sizes must be non-empty and sorted")
  sidx <- match(size_name, individual$ode_names)
  if (is.na(sidx))
    stop("size_name must be one of the ODE state names: ",
         paste(individual$ode_names, collapse = ", "))
  if (is.null(traits)) traits <- ff16_default_traits()

  # Pass 1: harvest the adaptive step schedule + initial state (the frozen schedule the
  # replay reproduces). grow_individual_bracket returns the step times and per-node states.
  brk <- grow_individual_bracket(individual, sizes, size_name, env, time_max, warn)
  y0  <- stats::setNames(individual$ode_state, individual$ode_names)
  pp  <- unlist(individual$strategy$pars)

  res <- ff16_grow_to_size_gradient_impl(pp, env, y0, brk$time, as.numeric(sizes),
                                         as.integer(sidx - 1L), traits, TRUE)
  rownames(res$time)  <- NULL
  res$sizes <- sizes
  res
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
