# TF24f stand-metric gradients (#472 scope B, build-order step 1). TF24f is the
# fast-acclimation TF24 variant whose optimal root-collar potential is a 6th ODE state
# tracked by gradient ascent (no per-step golden-section optimiser). Per
# notes/tf24-stand-gradient-scope.md, TF24f -- not TF24 -- is the right target for the
# stand CENSUS gradients: its leaf evaluation at the tracked collar is analytic / IFT-able,
# so the census number density's growth-rate-gradient term needs no curvature harvest. The
# R0 GATE below is the first deliverable: a double-precision census reconstruction that
# proves the collar-state replay is faithful (the prerequisite §7 flags) before the
# reverse-mode tape (R1) is built.

# Harvest a run-with-cache TF24f SCM into the pieces the census replay consumes (the
# TF24f mirror of tf24_harvest). Unlike TF24, no per-stage leaf-opt harvest is needed --
# the replay evaluates the analytic leaf at the tracked collar live -- so this returns
# only the frozen schedule + env, the cohort birth steps, the parameter vector and the
# acclimation knobs (k_acclim / use_ad_gradient). Requires shading_model = "crown-centre".
tf24f_harvest <- function(scm, species = 1L, birth_rate = NULL) {
  types <- extract_RcppR6_template_types(scm$parameters, "Parameters")
  if (!identical(types[[1]], "TF24f")) {
    stop("TF24f census gradients are implemented for the TF24f strategy only")
  }
  if (species < 1L || species > length(scm$patch$species)) {
    stop("species index out of range: stand has ", length(scm$patch$species),
         " species")
  }
  patch <- scm$patch
  sh <- patch$step_history
  eh <- patch$environment_history
  if (length(eh) < 1L) {
    stop("No resident schedule cached: run the SCM with control(save_RK45_cache = TRUE)")
  }
  sp    <- patch$species[[species]]
  nt    <- sp$node_times
  strat <- scm$parameters$strategies[[species]]
  pp    <- unlist(strat$pars)

  if (is.null(birth_rate)) {
    birth_rate <- scm$offspring_production[[species]] / scm$net_reproduction_ratios[[species]]
  }
  # Cohort birth steps (introductions land on step times).
  birth_step <- vapply(nt, function(t) which.min(abs(sh - t)) - 1L, integer(1))

  list(pp = pp, eh = eh, sh = sh, birth_step = birth_step, birth_rate = birth_rate,
       k_acclim = strat$k_acclim, use_ad_gradient = strat$use_ad_gradient, nt = nt)
}

# R0 gate (internal): double-precision census reconstruction of a TF24f resident stand.
# Re-evolves every cohort's {5 demog, tracked collar, log_density} over the frozen
# schedule and returns the reconstructed per-cohort heights / collar / log-densities and
# the census metric values. Used to confirm the collar-state replay reproduces the SCM's
# stored stand (heights / log_densities / opt_root_psi_state) and that the LAI reduction
# matches compute_competition(0), before the reverse-mode tape (R1) is added. Not yet a
# public gradient entry point -- it returns the recon, not d(metric)/d(theta).
tf24f_census_recon <- function(scm, metrics = c("LAI", "biomass", "size_moment"),
                               species = 1L, birth_rate = NULL) {
  h <- tf24f_harvest(scm, species, birth_rate)
  # The reconstruction's growth-rate gradient g' must match whatever the resident
  # SCM used (Node::growth_rate_gradient): the exact-AD path when the run set
  # control(node_gradient_exact_ad = TRUE), else the backward finite difference.
  exact_ad <- isTRUE(scm$parameters$strategies[[species]]$control$node_gradient_exact_ad)
  tf24f_census_recon_impl(h$pp, h$eh, h$sh, h$birth_step, h$birth_rate, h$k_acclim,
                          h$use_ad_gradient, metrics, exact_ad)
}

# TF24f frozen census trait gradient (#472 scope B, build-order step 2 -- R1 GATE).
# d(census metric)/d(theta) for the FROZEN (rare-mutant / invasion) resident light, by a
# central finite difference over the R0 census reconstruction: for each trait, perturb the
# parameter vector by +/- a relative step, re-run the double-precision collar-state replay
# against the SAME frozen resident environment, and difference the metric. The resident
# schedule + per-RK-stage light are held fixed (the invasion gradient), so this is the
# faithful finite-difference reference the reverse-mode AD tape (the actual R1) must
# reproduce -- and a usable prototype gradient in its own right. Returns a metrics x traits
# Jacobian and the reconstructed metric values. (The AD tape will replace the per-trait
# replays with one reverse sweep; this gate fixes the target it must hit.)
tf24f_census_gradient_fd <- function(scm, metrics = c("LAI", "biomass", "size_moment"),
                                     traits = NULL, species = 1L, birth_rate = NULL,
                                     rel_step = 1e-5) {
  if (is.null(traits)) traits <- tf24_default_traits()
  h <- tf24f_harvest(scm, species, birth_rate)
  exact_ad <- isTRUE(scm$parameters$strategies[[species]]$control$node_gradient_exact_ad)
  recon <- function(pp) tf24f_census_recon_impl(pp, h$eh, h$sh, h$birth_step,
    h$birth_rate, h$k_acclim, h$use_ad_gradient, metrics, exact_ad)$values

  values <- recon(h$pp)
  jac <- matrix(0, length(metrics), length(traits),
                dimnames = list(metrics, traits))
  for (j in seq_along(traits)) {
    tr <- traits[j]
    if (!tr %in% names(h$pp)) stop("unknown TF24f trait: ", tr)
    d <- rel_step * max(abs(h$pp[[tr]]), 1e-8)
    pp_p <- h$pp; pp_p[[tr]] <- pp_p[[tr]] + d
    pp_m <- h$pp; pp_m[[tr]] <- pp_m[[tr]] - d
    jac[, j] <- (recon(pp_p) - recon(pp_m)) / (2 * d)
  }
  list(jacobian = jac, values = values)
}

# TF24f individual grow-to-size trait gradient (#472 scope B, the "individuals" surface
# -- prototype). d(t*)/d(theta) and d(state at t*)/d(theta) for a single plant grown in a
# FIXED environment to target size(s), by central finite difference over
# grow_individual_to_size. No resident feedback (the env is given) and no canopy/density,
# so this is the lightest gradient surface; for TF24f the tracked collar is re-evolved
# inside each grow (it is one of the ODE states), so the FD captures the collar's response
# automatically -- which is exactly why an exact AD version is the heavier follow-up (the
# tracked collar is strongly theta-dependent; see notes/tf24-stand-gradient-scope.md). The
# FD here is the prototype + the reference that AD version must reproduce. Traits are
# perturbed on the (post-hyperpar) strategy parameters directly, matching the census FD.
tf24f_grow_individual_to_size_gradient_fd <- function(individual, sizes, size_name, env,
                                                      traits = NULL, time_max = Inf,
                                                      warn = FALSE, rel_step = 1e-5) {
  if (!grepl("^TF24f", individual$strategy_name))
    stop("tf24f_grow_individual_to_size_gradient_fd is for the TF24f strategy only")
  if (is.null(traits)) traits <- tf24_default_traits()

  grow <- function(ind) {
    r <- grow_individual_to_size(ind, sizes, size_name, env, time_max, warn)
    list(time = r$time, state = r$state)
  }
  perturb <- function(tr, delta) {
    s <- individual$strategy
    pars <- s$pars
    if (!tr %in% names(pars)) stop("unknown TF24f trait: ", tr)
    pars[[tr]] <- pars[[tr]] + delta
    s$pars <- pars
    grow(TF24f_Individual(s))
  }

  base <- grow(individual)
  nS <- length(sizes); comp <- colnames(base$state); nC <- length(comp)
  d_time <- matrix(0, nS, length(traits), dimnames = list(NULL, traits))
  d_state <- array(0, c(nS, nC, length(traits)), dimnames = list(NULL, comp, traits))
  s0 <- individual$strategy$pars
  for (j in seq_along(traits)) {
    tr <- traits[j]
    d <- rel_step * max(abs(s0[[tr]]), 1e-8)
    rp <- perturb(tr, d); rm <- perturb(tr, -d)
    d_time[, j] <- (rp$time - rm$time) / (2 * d)
    d_state[, , j] <- (rp$state - rm$state) / (2 * d)
  }
  list(sizes = sizes, time = base$time, state = base$state,
       d_time = d_time, d_state = d_state)
}

# TF24f RESIDENT (coupled) census trait gradient (#472 scope B, the resident-feedback
# surface -- prototype). The TOTAL stand-level d(census metric)/d(theta): unlike the
# frozen (rare-mutant) gradient, every cohort's height + density feeds back through the
# canopy light that the whole stand reads, so the feedback routinely dominates and can
# flip the sign relative to the frozen reading. Computed as a central finite difference
# over the FULL SCM: perturb a (post-hyperpar) strategy parameter, re-run run_scm on the
# SAME fixed node schedule, and difference the realised stand metric. This is the ground-
# truth resident gradient (the SCM responds in full) and the reference the coupled AD
# engine must reproduce; it is slow (one SCM solve per trait per side), so keep `traits`
# small. LAI is read from the realised patch (compute_competition(0)); size_moment is the
# size-distribution first moment Sum density_i * height_i (trapezium over the stand).
tf24f_resident_census_gradient_fd <- function(p, env, ctrl,
                                              metrics = c("LAI", "size_moment"),
                                              traits = NULL, species = 1L,
                                              rel_step = 1e-4) {
  if (is.null(traits)) traits <- tf24_default_traits()
  metric_set <- c("LAI", "size_moment")
  bad <- setdiff(metrics, metric_set)
  if (length(bad)) stop("unsupported resident metric(s): ", paste(bad, collapse = ", "))

  # Realised stand metrics from a completed SCM (the coupled, self-shaded stand).
  stand_metrics <- function(scm) {
    patch <- scm$patch
    out <- c(LAI = patch$compute_competition(0))
    sp <- patch$species[[species]]
    h <- sp$heights; dens <- exp(sp$log_densities)
    ord <- order(h, decreasing = TRUE)
    h <- h[ord]; dens <- dens[ord]
    # size_moment = trapezium of density*height over descending heights down to h0.
    phi <- dens * h
    sm <- 0
    if (length(h) > 1) sm <- sum(0.5 * (h[-length(h)] - h[-1]) * (phi[-length(h)] + phi[-1]))
    out["size_moment"] <- sm
    out[metrics]
  }
  run_with <- function(pp_override) {
    pmod <- p
    if (!is.null(pp_override)) {
      s <- pmod$strategies[[species]]
      s$pars <- pp_override
      pmod$strategies[[species]] <- s
    }
    run_scm(pmod, env, ctrl, refine_schedule = FALSE)
  }

  base_scm <- run_with(NULL)
  values <- stand_metrics(base_scm)
  pars0 <- p$strategies[[species]]$pars
  jac <- matrix(0, length(metrics), length(traits),
                dimnames = list(metrics, traits))
  for (j in seq_along(traits)) {
    tr <- traits[j]
    if (!tr %in% names(pars0)) stop("unknown TF24f trait: ", tr)
    d <- rel_step * max(abs(pars0[[tr]]), 1e-8)
    pp_p <- pars0; pp_p[[tr]] <- pp_p[[tr]] + d
    pp_m <- pars0; pp_m[[tr]] <- pp_m[[tr]] - d
    mp <- stand_metrics(run_with(pp_p)); mm <- stand_metrics(run_with(pp_m))
    jac[, j] <- (mp - mm) / (2 * d)
  }
  list(jacobian = jac, values = values)
}
