##' Construct a \code{Control} object. \code{control()} is a lowercase alias for
##' the \code{Control()} constructor, whose defaults are the pragmatic,
##' fast-ish settings used for essentially all of plant's runs (see
##' \code{control.cpp}). \code{control_accurate()} tightens the ODE and schedule
##' tolerances for high-accuracy runs at the cost of speed. \code{control_tf24()}
##' sets TF24's step control: the relative tolerance \code{tol} with the absolute
##' one at 1e-4 of it, the soil layers' error weight at 10, each weight bounded
##' at 100, steps capped at 15 days, and the soil stepped alone where the stand
##' draws under a tenth of the water moving through it. The window's weights are
##' \code{control_window()}'s.
##'
##' \code{control_window()} sets the window's weights from a pilot of the
##' analysis: a step starting at \code{t} takes the factor
##' \code{1 / min(max(R(t) / R0, r_min), 1)}, where \code{R(t)} is the largest
##' share of offspring production still to be earned after \code{t}, over the
##' pilot's stand and each invader walked on its recording. A pilot of 54 uniform
##' introductions at a tolerance of 1e-3 reads \code{R} within 10\% wherever it
##' is at least 1e-3. Under constant rain a uniform pilot lumps the founders, and
##' a grid that resolves them reads it within 1.1\%.
##'
##' The SCM's adaptive ODE stepper multiplies each state's error level by a
##' weight, which decides which steps it takes and nothing else:
##' \code{ode_weight_soil} on the soil layers, \code{ode_weight_accumulator} on
##' the flux accumulators, and on every state, for a step starting at time
##' \code{t}, the entry of \code{ode_weight_factors} paired with the last of
##' \code{ode_weight_times} at or before \code{t}. The times start at 0 and are
##' sorted. \code{ode_weight_max} bounds each state's weight once the factor
##' multiplies it. The weights default to 1, the schedule to empty and the bound
##' to \code{Inf}, which leaves every step as it was.
##'
##' \code{ode_split_sign_changes} integrates each TF24 node in pieces between the
##' sign changes of its net production inside a step, in the field sampled at five
##' fractions of the step. It is off by default, which leaves every run as it was.
##'
##' Under \code{control_tf24()} the bound keeps the soil's own error test in force:
##' unbounded, the soil's weight times a window factor of 100 reaches 1000, an
##' accepted step can carry a soil stage to the potential ceiling, and TF24's
##' gradient is refused there. The cap keeps an invader's storage pools stable,
##' which Cash-Karp loses on steps past 26 days.
##'
##' @title Control presets
##' @param ... Named control fields, passed to \code{Control()}.
##' @param base An optional \code{Control} object to start from; defaults are
##'   used if omitted.
##' @return A \code{Control} object.
##' @rdname control_presets
##' @export
control <- function(...) Control(...)

##' @rdname control_presets
##' @export
control_accurate <- function(base = Control()) {
  base$ode_tol_rel       <- 1e-6
  base$ode_tol_abs       <- 1e-6
  base$ode_step_size_max <- 1e-1
  base$schedule_eps      <- 1e-3
  base
}

##' @rdname control_presets
##' @param tol The relative tolerance, \code{ode_tol_rel}.
##' @export
control_tf24 <- function(tol = 3e-5, base = Control()) {
  base$ode_tol_rel       <- tol
  base$ode_tol_abs       <- 1e-4 * tol
  base$ode_weight_soil   <- 10
  base$ode_weight_max    <- 100
  base$ode_step_size_max <- 15 / 365
  base$ode_soil_alone_share <- 0.1
  base
}

##' @rdname control_presets
##' @param pilot A run of the analysis on the birth-date coordinate, of one
##'   species. Each invader is walked on its recording, so it ends holding the
##'   last walk.
##' @param invaders A list of \code{Parameters}, each walked on the pilot's
##'   introductions.
##' @param R0 The share still to be earned below which a step's weight rises
##'   from 1.
##' @param r_min The share of \code{R0} at which the weight stops rising, at
##'   \code{1 / r_min}.
##' @export
control_window <- function(pilot, invaders = list(), base = Control(),
                           R0 = 0.1, r_min = 0.01) {
  species <- pilot$patch$species
  if (length(species) != 1 || !species[[1]]$density_in_birth_date) {
    stop("control_window() reads a pilot of one species on the birth-date coordinate")
  }
  ## The share still to be earned after each row the last run or walk recorded,
  ## each node weighted as offspring production weights it.
  to_earn <- function() {
    rows <- pilot$store_trajectory()
    sp <- pilot$patch$species[[1]]
    per <- sp$ode_size / sp$size
    k <- match("offspring_produced_survival_weighted", sp$new_node$ode_names)
    environment_size <- length(rows[[length(rows)]]$state) - sp$ode_size
    weight <- head(sp$establishment_weights, -1) * sp$patch_densities *
      sp$extrinsic_drivers$evaluate_range("birth_rate", sp$node_times)
    earned <- vapply(rows, function(row) {
      j <- seq_len((length(row$state) - environment_size) / per)
      sum(weight[j] * row$state[per * (j - 1) + k])
    }, 0)
    list(time = vapply(rows, `[[`, 0, "time"),
         R = 1 - earned / earned[length(earned)])
  }
  times <- pilot$parameters$node_schedule_times
  stand <- to_earn()
  R <- stand$R
  for (p in invaders) {
    p$node_schedule_times <- times
    pilot$run_mutant(p)
    walk <- to_earn()
    ## An invader that earns nothing has no share to protect.
    R <- pmax(R, walk$R[findInterval(stand$time, walk$time)], na.rm = TRUE)
  }
  base$ode_weight_times <- stand$time
  base$ode_weight_factors <- 1 / pmin(pmax(R / R0, r_min), 1)
  base
}


##' Basic default settings for a given strategy, environment only really
##' used for templating initially and will be overloaded later by passing
##' an environment to the SCM API (suggesting perhaps the template could be
##' removed).
##' @title Basic default parameters for a given strategy
##' @author Rich FitzJohn
##' @param type Any strategy name as a string, e.g.: \code{"FF16"}.
##' @param env And environment object
##' @export
scm_base_parameters <- function(type = NA, env = environment_type(type)) {
  Parameters(type, env)()
}


##' Run the SCM.
##'
##' The node-introduction schedule can be adaptively refined in C++ by setting
##' \code{refine_schedule = TRUE} (this replaces the former \code{build_schedule}
##' function). Setting \code{collect = TRUE} returns tidied output collected at
##' every ODE step (replacing the former \code{run_scm_collect}); otherwise the
##' \code{SCM} object itself is returned for interrogation.
##'
##' @title Run SCM
##' @param p Parameters object
##' @param env Environment object (defaults to the strategy's environment)
##' @param ctrl Control object
##' @param events An \code{\link{events}} object giving the discrete events to
##'   apply during the run — rainfall pulses, harvest, and node introductions
##'   themselves. When \code{NULL} (the default) the schedule is taken from
##'   \code{p$node_schedule_times}, i.e. introductions only.
##' @param refine_schedule Should the node-introduction schedule be adaptively
##'   refined before/while running (using \code{schedule_eps} and
##'   \code{schedule_nsteps} from \code{ctrl})? Refinement records the ODE
##'   schedule its final run took into \code{p$ode_times},
##'   \code{p$ode_step_sizes}, \code{p$ode_alone_slopes} and
##'   \code{p$ode_alone_steps}, so a later run of those parameters replays it
##'   exactly rather than choosing its own steps again.
##' @param record_trajectory Should the run keep the state at every accepted
##'   step? A gradient sweeps those states and cannot recover them from a
##'   finished run, so a run that did not keep them is repeated -- one whole
##'   forward integration. Asked for here because the flag has to be set before
##'   the run that fills it, and this function is where that run happens. With
##'   \code{refine_schedule = TRUE} the refinement's own runs do not keep them;
##'   only the final run does.
##' @param collect Should tidied results be collected at every step and
##'   returned (instead of the \code{SCM} object)?
##' @return When \code{collect = FALSE}, an \code{SCM} object. When
##'   \code{collect = TRUE}, a list of tidied patch output with
##'   \code{offspring_production}, \code{net_reproduction_ratios} and the
##'   (possibly refined) parameters \code{p}.
##' @author Rich FitzJohn
##' @rdname run_scm
##' @export
run_scm <- function(p, env = NULL,
                    ctrl = control(),
                    refine_schedule = FALSE, collect = FALSE,
                    record_trajectory = FALSE, events = NULL) {

  types <- extract_RcppR6_template_types(p, "Parameters")

  if (is.null(env))
    env <- Environment(types[[1]])

  # An ODE schedule carried by the parameters is taken, because a schedule is
  # there to be used: p$ode_times with p$ode_step_sizes and the soil's records
  # replays a recorded run exactly, and p$ode_times alone stops at a grid the
  # caller chose. To integrate freely, carry neither.

  ## No events supplied: the schedule comes from p$node_schedule_times, as it
  ## did before events existed. An empty Events object is how that is signalled
  ## across the boundary.
  if (is.null(events))
    events <- empty_events()

  scm <- do.call('SCM', types)(p, env, events, ctrl)
  if (collect) {
    scm$collect <- TRUE
  }

  if (refine_schedule) {
    # Refinement's runs are bisected against and discarded, so they keep no
    # states; the run a sweep walks is one more, after the schedule settles.
    scm$refine_schedule()
    if (record_trajectory) {
      scm$record_trajectory <- TRUE
      scm$run()
    }
  } else {
    scm$record_trajectory <- record_trajectory
    scm$run()
  }

  if (!collect) {
    return(scm)
  }

  results <- lapply(scm$history, "[[", "state") |> tidy_patch()
  results[["offspring_production"]] <- scm$offspring_production
  results[["net_reproduction_ratios"]] <- scm$net_reproduction_ratios
  results[["p"]] <- scm$parameters
  ## Events are supplied separately from `p`, so `p` alone does not describe the
  ## run. Carry both the requested schedule and what was actually applied, or a
  ## collected result silently loses the whole event record -- including how much
  ## of each pulse the soil accepted and how much it shed.
  results[["events"]] <- scm$events
  results[["event_log"]] <- scm$event_log

  results
}

##' Export the full state of a patch from a (run) \code{SCM} so it can be
##' re-imported to seed a new run (see \code{\link{set_initial_state}}). The
##' exported object captures everything needed to reproduce the patch's forward
##' trajectory: every node's ODE state plus the per-node birth bookkeeping
##' (introduction time, patch-age density and survival probability at birth)
##' that is not part of the ODE state but feeds the rates and lifetime-fitness
##' integrals, the patch age, and the not-yet-introduced ("residual") portion of
##' the node-introduction schedule.
##'
##' @title Export patch state from an SCM
##' @param scm An \code{SCM} object that has been run.
##' @param step Optional 1-based index into \code{scm$history} (requires the run
##'   to have been performed with \code{collect = TRUE}); when \code{NULL}
##'   (default) the SCM's current/final patch is used.
##' @return A list describing the patch state, suitable for
##'   \code{\link{set_initial_state}}: \code{time}, \code{n} (nodes per species),
##'   \code{ode_state} (flat), the per-species lists \code{node_times},
##'   \code{patch_density} and \code{pr_patch_survival}, and the residual
##'   \code{node_schedule_times}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{run_scm}}
##' @export
export_patch_state <- function(scm, step = NULL) {
  patch <- if (is.null(step)) scm$patch else scm$history[[step]]
  time <- patch$time
  species <- patch$species

  ## Residual schedule: the introductions not yet represented among the seeded
  ## nodes. A patch snapshot sits integrated *to* its next introduction time
  ## without having introduced that node yet (its node times are all strictly
  ## below `time`), so the residual keeps every original time at or after `time`
  ## -- the resumed run introduces the node due at `time` itself, then the rest.
  tol <- 1e-8
  residual <- lapply(scm$parameters$node_schedule_times,
                     function(tt) tt[tt >= time - tol])

  list(
    time = time,
    n = as.integer(vapply(species, function(s) s$size, numeric(1))),
    ode_state = patch$ode_state,
    node_times = lapply(species, function(s) s$node_times),
    patch_density = lapply(species, function(s) s$patch_densities),
    pr_patch_survival = lapply(species, function(s) s$pr_patch_survival_at_birth),
    node_schedule_times = residual
  )
}

##' Build an initial patch state from a specified size distribution, for seeding
##' a patch at age 0 with pre-existing plants instead of growing it from empty
##' (the ecological motivation of \code{plant} issue #304). The returned object
##' is consumed by \code{\link{set_initial_state}}.
##'
##' Each species' initial nodes are described by their heights and densities;
##' the remaining ODE state (mortality, fecundity, accumulated reproduction,
##' heartwood, ...) starts at zero, all nodes are introduced at patch age 0, and
##' the recruitment schedule continues for \code{t > 0} (the seeded distribution
##' replaces the \code{t = 0} recruit). Pathologically large/dense initial
##' conditions can produce non-finite densities; \code{run_scm} guards against
##' this and errors with a suggestion to use more plausible inputs.
##'
##' @title Build an initial size distribution
##' @param p A \code{Parameters} object.
##' @param heights Per-species node heights: a numeric vector (single species) or
##'   a list of numeric vectors (one per species).
##' @param densities Per-species node densities (same shape as \code{heights}).
##'   Supply this or \code{log_densities}.
##' @param log_densities Per-species node log-densities (alternative to
##'   \code{densities}).
##' @param env Environment object (defaults to the strategy's environment).
##' @param ctrl Control object.
##' @param birth_dates Per-species node birth dates, before 0 and ascending as
##'   the heights descend (same shape as \code{heights}). Required on the
##'   birth-date coordinate, where each node's density is its birth rate times
##'   its survival and \code{densities} are densities in birth date: the
##'   survival is set to give that density, and each interval's establishment to
##'   its width, as if every seed established. Elsewhere every node is born at 0.
##' @return A state list suitable for \code{\link{set_initial_state}}.
##' @seealso \code{\link{set_initial_state}}, \code{\link{export_patch_state}}
##' @export
make_initial_state <- function(p, heights, densities = NULL,
                               log_densities = NULL, env = NULL,
                               ctrl = control(), birth_dates = NULL) {
  types <- extract_RcppR6_template_types(p, "Parameters")
  if (is.null(env)) {
    env <- Environment(types[[1]])
  }
  n_spp <- length(p$strategies)

  as_list <- function(z) if (is.list(z)) z else list(z)
  heights <- as_list(heights)
  if (length(heights) != n_spp) {
    stop("`heights` must have one entry per species (", n_spp, ")")
  }
  if (is.null(log_densities)) {
    if (is.null(densities)) {
      stop("supply either `densities` or `log_densities`")
    }
    log_densities <- lapply(as_list(densities), log)
  } else {
    log_densities <- as_list(log_densities)
  }

  ## Learn the node ODE layout and the (node-free) environment ODE tail from a
  ## fresh patch, plus the patch-age-0 disturbance weights for birth bookkeeping.
  patch <- do.call("Patch", types)(p, env, ctrl)
  ode_names <- patch$species[[1]]$new_node$ode_names
  hi <- match("height", ode_names)
  birth_date <- ctrl$node_density_in_birth_date
  ldi <- match(if (birth_date) "mortality" else "log_density", ode_names)
  ii <- match(c("interval_establishment", "interval_establishment_moment"),
              ode_names)
  if (is.na(hi) || is.na(ldi) || (birth_date && anyNA(ii))) {
    stop("could not locate the node's height and density in its ODE names")
  }
  if (birth_date) {
    if (is.null(birth_dates)) {
      stop("the birth-date coordinate needs each node's `birth_dates`")
    }
    birth_dates <- as_list(birth_dates)
  }
  node_ode_size <- length(ode_names)
  env_state <- patch$ode_state # fresh patch has no nodes: environment ODE only
  pr_surv0 <- patch$pr_survival(0)
  dens0 <- patch$density(0)

  n <- integer(n_spp)
  ode_chunks <- vector("list", n_spp)
  node_times <- vector("list", n_spp)
  patch_density <- vector("list", n_spp)
  pr_patch_survival <- vector("list", n_spp)
  for (i in seq_len(n_spp)) {
    h <- heights[[i]]
    ld <- log_densities[[i]]
    if (length(h) != length(ld)) {
      stop("heights and densities must have equal length (species ", i, ")")
    }
    ## The model requires nodes ordered by decreasing height.
    o <- order(h, decreasing = TRUE)
    h <- h[o]
    ld <- ld[o]
    mat <- matrix(0, nrow = node_ode_size, ncol = length(h))
    mat[hi, ] <- h
    times <- rep(0, length(h))
    if (birth_date) {
      times <- birth_dates[[i]][o]
      width <- diff(c(times, 0))
      if (length(times) != length(h) || any(width <= 0)) {
        stop("birth_dates must be before 0 and ascend as the heights descend ",
             "(species ", i, ")")
      }
      drivers <- patch$species[[i]]$extrinsic_drivers
      birth_rate <- vapply(times, function(t) drivers$evaluate("birth_rate", t),
                           numeric(1))
      mat[ldi, ] <- log(birth_rate) - ld
      mat[ii, ] <- rbind(width, width^2 / 2)
    } else {
      mat[ldi, ] <- ld
    }
    ode_chunks[[i]] <- as.vector(mat) # column-major: node-by-node, matching set_ode_state
    n[i] <- length(h)
    node_times[[i]] <- times
    patch_density[[i]] <- rep(dens0, length(h))
    pr_patch_survival[[i]] <- rep(pr_surv0, length(h))
  }

  list(
    time = 0,
    n = n,
    ode_state = c(unlist(ode_chunks, use.names = FALSE), env_state),
    node_times = node_times,
    patch_density = patch_density,
    pr_patch_survival = pr_patch_survival,
    ## Continue recruitment for t > 0; the seeded distribution replaces the t=0 recruit.
    node_schedule_times = lapply(p$node_schedule_times, function(tt) tt[tt > 1e-8])
  )
}

##' Write an exported patch state (from \code{\link{export_patch_state}}) into a
##' \code{Parameters} object so that the next \code{\link{run_scm}} starts from
##' that state instead of an empty patch. The state is carried on the
##' \code{Parameters} object (rather than passed separately) so the run stays
##' self-describing and reproducible, and so the seeding survives the reset at
##' the start of every run / schedule refinement.
##'
##' @title Seed Parameters with an initial patch state
##' @param p A \code{Parameters} object (its strategies must match the exported
##'   state).
##' @param state An exported state list from \code{\link{export_patch_state}}.
##' @return The modified \code{Parameters} object.
##' @seealso \code{\link{export_patch_state}}, \code{\link{run_scm}}
##' @export
set_initial_state <- function(p, state) {
  n_spp <- length(p$strategies)
  if (length(state$n) != n_spp) {
    stop("State has ", length(state$n),
         " species but Parameters has ", n_spp, " strategies")
  }
  p$initial_state <- state$ode_state
  p$n_initial_cohorts <- as.integer(state$n)
  p$initial_node_times <- unlist(state$node_times, use.names = FALSE)
  p$initial_patch_density <- unlist(state$patch_density, use.names = FALSE)
  p$initial_pr_patch_survival <- unlist(state$pr_patch_survival, use.names = FALSE)
  p$initial_time <- state$time
  if (!is.null(state$node_schedule_times)) {
    p$node_schedule_times <- state$node_schedule_times
  }
  p
}

