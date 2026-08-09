## An independent finite-difference reference for one cohort's rate vector.
##
## A reverse-mode sensitivity pass over a stand has to transpose one cohort's
## physiology at one solver stage. A transpose checked only against the forward
## code it was derived from passes whenever both share an error, so the check
## needs a reference built by a different route. A cohort's own state is six
## numbers, which is few enough that the whole Jacobian can be formed
## exhaustively by central finite difference -- one perturbation per state
## component -- at a cost of 12 rate evaluations.
##
## Nothing here is pinned. The reference is recomputed on demand: this family's
## numbers are not bit-reproducible across platforms, so a recorded table would
## detect change rather than incorrectness. The one part that IS asserted
## exactly is the sparsity pattern (see `cohort_jacobian_pattern`), which
## follows from which states the rate code reads and is therefore identically
## zero rather than small.

##' Rate vector of a single cohort at one solver stage.
##'
##' Everything a `Node` produces from its state in one right-hand-side
##' evaluation: the strategy's ODE rates, the two node-level rates
##' (survival-weighted offspring production and the density rate), and the
##' per-resource consumption rates that feed the patch resource balance. For
##' TF24 the last of these is per-soil-layer water uptake in m/yr, per unit
##' density.
##'
##' This reads the node at whatever state it is already in; it recomputes the
##' node's rates as a side effect but does not move its state.
##'
##' @param node A `Node` object, e.g. from `patch$species[[i]]$node_at(j)`.
##' @param environment The environment to evaluate in, held fixed. This is a
##'   partial derivative: perturbing a cohort in a real patch would move the
##'   canopy and the soil, and that feedback is deliberately excluded.
##' @param pr_patch_survival Probability of patch survival to the current time,
##'   i.e. `patch$pr_survival(patch$time)`.
##' @return A named numeric vector.
##' @export
cohort_rate_vector <- function(node, environment, pr_patch_survival) {
  node$compute_rates(environment, pr_patch_survival)
  internals <- node$individual$internals
  out <- c(node$ode_rates, internals$consumption_rates)
  names(out) <- cohort_rate_names(node)
  out
}

##' @rdname cohort_rate_vector
##' @export
cohort_rate_names <- function(node) {
  n_resources <- node$individual$internals$resource_size
  c(node$ode_names,
    if (n_resources > 0) paste0("consumption_", seq_len(n_resources)))
}

##' Jacobian of a cohort's rate vector with respect to its own state.
##'
##' Formed by central finite difference, one perturbation per state component.
##' The columns are the strategy's ODE states (for TF24: height, mortality,
##' fecundity, the two heartwood accumulators and the storage reserve); the
##' rows are `cohort_rate_names(node)`. The two node-level accumulators
##' (offspring produced, log density) are not perturbed: neither is read by any
##' rate, and the density scaling they carry is a property of the node rather
##' than of the physiology being transposed.
##'
##' The step is relative: component `i` is moved by `step * abs(x[i])`, or by
##' `step` when `x[i]` is zero. Do not pick a step by eye. The default is the
##' geometric centre of the range over which every rate that depends on the
##' physiology holds still, measured with `cohort_jacobian_sweep` and
##' `cohort_jacobian_plateau` at five nodes of a two-year TF24 stand: the range
##' common to all of them is 1.8e-7 to 5.6e-5 at a flatness of 1e-6. Its lower
##' edge is not machine precision but the leaf sub-model's own solver
##' tolerances, which leave a relative noise floor of 5e-10 to 4e-9 on the
##' rates; its upper edge is curvature, where the truncation error starts
##' growing as the square of the step.
##'
##' The density rate is the exception and does not share that plateau: see
##' `cohort_jacobian_plateau`.
##'
##' Perturbation goes through `Node$ode_state`, which routes to
##' `Individual::set_state` and so refreshes the aux slots derived from height
##' (`competition_effect`, `height_inverse`). `Internals$set_state` is also
##' exposed to R and does not: it writes the state and leaves those two holding
##' the previous height's values, which then feed the whole rate path. Perturb
##' through the node or the individual, not through the internals.
##'
##' @inheritParams cohort_rate_vector
##' @param step Relative finite-difference step.
##' @return A matrix with `length(cohort_rate_names(node))` rows and
##'   `node$individual$ode_size` columns, with dimnames.
##' @export
cohort_rate_jacobian <- function(node, environment, pr_patch_survival,
                                 step = 3e-6) {
  x0 <- node$ode_state
  n_state <- node$individual$ode_size
  state_names <- node$individual$ode_names
  on.exit({
    node$ode_state <- x0
    node$compute_rates(environment, pr_patch_survival)
  }, add = TRUE)

  f <- function(x) {
    node$ode_state <- x
    cohort_rate_vector(node, environment, pr_patch_survival)
  }

  columns <- lapply(seq_len(n_state), function(i) {
    dx <- step * (if (x0[[i]] == 0) 1 else abs(x0[[i]]))
    up <- x0; up[[i]] <- x0[[i]] + dx
    down <- x0; down[[i]] <- x0[[i]] - dx
    (f(up) - f(down)) / (2 * dx)
  })

  jacobian <- do.call(cbind, columns)
  dimnames(jacobian) <- list(cohort_rate_names(node), state_names)
  jacobian
}

##' Sweep the finite-difference step over a grid.
##'
##' Used to locate the plateau on which the difference quotient is a derivative
##' rather than either truncation error or subtractive cancellation.
##'
##' @inheritParams cohort_rate_jacobian
##' @param steps Relative steps to sweep, in increasing order.
##' @return A 3-d array indexed `[rate, state, step]`, with `steps` recorded as
##'   the third dimnames component and as the `steps` attribute.
##' @export
cohort_jacobian_sweep <- function(node, environment, pr_patch_survival,
                                  steps = 10^seq(-11, -2, by = 0.25)) {
  jacobians <- lapply(steps, function(s) {
    cohort_rate_jacobian(node, environment, pr_patch_survival, step = s)
  })
  out <- array(unlist(jacobians),
               dim = c(dim(jacobians[[1]]), length(steps)),
               dimnames = c(dimnames(jacobians[[1]]),
                            list(format(steps, digits = 3))))
  attr(out, "steps") <- steps
  out
}

##' Plateau edges of a finite-difference step sweep.
##'
##' For each entry of the Jacobian, the longest run of consecutive steps over
##' which the difference quotient changes by less than `tol` in relative terms.
##' Both edges are reported: the lower one is where cancellation against the
##' solver's own tolerances takes over, the upper one where curvature does.
##' Entries whose magnitude never exceeds `zero_tol` are skipped -- a
##' structurally zero entry has no plateau to find.
##'
##' The `log_density` rows do not behave like the rest and should be read
##' separately. `Node::compute_rates` builds the density rate from a backward
##' difference of the height growth rate taken at a fixed absolute step
##' (`control()$node_gradient_eps`, 1e-6), so the density rate carries a
##' relative noise floor of 4e-8 where the other rates carry 5e-10 to 4e-9.
##' At a flatness of 1e-6 it has no plateau at all; at 1e-4 it has one, but a
##' decade or two coarser than everything else's.
##'
##' @param sweep Output of `cohort_jacobian_sweep`.
##' @param tol Relative change between adjacent steps that still counts as flat.
##' @param zero_tol Entries never exceeding this magnitude are skipped.
##' @return A data frame with one row per surviving entry: `rate`, `state`, the
##'   plateau's `lower` and `upper` step edges, its width in decades, the
##'   geometric-centre step `centre`, and the `value` there.
##' @export
cohort_jacobian_plateau <- function(sweep, tol = 1e-6, zero_tol = 0) {
  steps <- attr(sweep, "steps")
  rates <- dimnames(sweep)[[1]]
  states <- dimnames(sweep)[[2]]

  rows <- list()
  for (i in seq_along(rates)) {
    for (j in seq_along(states)) {
      v <- sweep[i, j, ]
      if (max(abs(v)) <= zero_tol) {
        next
      }
      run <- longest_flat_run(v, tol)
      if (is.null(run)) {
        next
      }
      centre <- exp(mean(log(steps[run])))
      k <- which.min(abs(log(steps[run]) - log(centre)))
      rows[[length(rows) + 1L]] <- data.frame(
        rate = rates[[i]], state = states[[j]],
        lower = steps[[run[[1]]]], upper = steps[[run[[length(run)]]]],
        decades = log10(steps[[run[[length(run)]]]] / steps[[run[[1]]]]),
        centre = steps[run][[k]], value = v[run][[k]],
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(rate = character(0), state = character(0),
                      lower = numeric(0), upper = numeric(0),
                      decades = numeric(0), centre = numeric(0),
                      value = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

## Longest run of consecutive entries of `v` whose successive relative change
## stays below `tol`. Returns the indices, or NULL if no two adjacent entries
## agree that closely.
longest_flat_run <- function(v, tol) {
  n <- length(v)
  if (n < 2L) {
    return(NULL)
  }
  scale <- pmax(abs(v[-n]), abs(v[-1L]))
  flat <- abs(v[-1L] - v[-n]) <= tol * scale
  flat[!is.finite(flat)] <- FALSE
  if (!any(flat)) {
    return(NULL)
  }
  r <- rle(flat)
  best <- which(r$values & r$lengths == max(r$lengths[r$values]))[[1]]
  end <- cumsum(r$lengths)[[best]]
  start <- end - r$lengths[[best]] + 1L
  seq.int(start, end + 1L)
}

##' Entries of the cohort Jacobian that are identically zero.
##'
##' The pattern is a statement about the model, not a measurement: it follows
##' from which state components the rate code reads at all. It is asserted with
##' exact comparison, not a tolerance, and the transpose of this pattern is what
##' any future adjoint has to satisfy.
##'
##' For TF24 the six states are read in exactly these places
##' (`src/tf24_strategy.cpp`, `inst/include/plant/node.h`):
##'
##' * `height` -- read throughout `compute_rates` and
##'   `net_mass_production_dt`. No zeros in this column.
##' * `mortality` -- reaches `mortality_dt` only through a
##'   `util::is_finite` test, so the mortality rate does not vary with it at
##'   all; the one rate that does read it is the node's survival-weighted
##'   offspring production, through `exp(-mortality)`.
##' * `fecundity` -- set as a rate and never read back. The whole column is
##'   zero.
##' * `area_heartwood` -- read only by `area_stem`/`diameter_stem`, which the
##'   rate path does not call. The whole column is zero.
##' * `mass_heartwood` -- read only by `mass_total`/`mass_above_ground`, same.
##'   The whole column is zero.
##' * `storage` -- read once, at the top of the storage block, and so reaches
##'   height, fecundity, storage and mortality. It does not reach the two
##'   heartwood rates, which are turnover-driven and depend on height alone;
##'   and it does not reach water consumption, because the leaf is solved by
##'   `net_mass_production_dt`, which is not passed the reserve.
##'
##' @param node A `Node` object, used for its rate and state names.
##' @return A logical matrix shaped like `cohort_rate_jacobian`, `TRUE` where
##'   the entry must be identically zero.
##' @export
cohort_jacobian_pattern <- function(node) {
  rates <- cohort_rate_names(node)
  states <- node$individual$ode_names
  if (node$individual$strategy_name != "TF24") {
    stop("cohort_jacobian_pattern: pattern derived for TF24 only, got ",
         node$individual$strategy_name)
  }

  zero <- matrix(FALSE, length(rates), length(states),
                 dimnames = list(rates, states))
  consumption <- grep("^consumption_", rates)

  zero[, "fecundity"] <- TRUE
  zero[, "area_heartwood"] <- TRUE
  zero[, "mass_heartwood"] <- TRUE

  zero[, "mortality"] <- TRUE
  zero["offspring_produced_survival_weighted", "mortality"] <- FALSE

  zero[c("area_heartwood", "mass_heartwood"), "storage"] <- TRUE
  zero[consumption, "storage"] <- TRUE

  zero
}
