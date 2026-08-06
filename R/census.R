## A stand census, and the first of the two terms in its sensitivity to a
## strategy parameter.
##
## A census is a quadrature over the size distribution: total leaf area, stem
## (basal) area and above-ground mass are each a sum over nodes of a quadrature
## weight times a stem number density times a per-individual quantity.
##
## Differentiating a census with respect to a parameter gives two terms. The
## parameter changes how the stand develops, so the size distribution moves; and
## the parameter changes what each individual *is* at whatever size it happens
## to have reached. Only the second -- the DIRECT term, taken at fixed state --
## is computed here. See inst/include/plant/census.h.

##' Names of the implemented stand census metrics.
##'
##' Each metric is a quadrature over the size distribution of a different
##' per-individual quantity.
##'
##' @return A character vector of metric names.
##' @export
census_metrics <- function() {
  census_metric_names()
}

##' Registered parameters for the stand census direct sensitivity term.
##'
##' The union, over the census metrics, of the strategy parameters the census
##' algebra reads. These are strategy parameters (`s$pars$...`), not traits: a
##' trait sensitivity needs the hyperpar Jacobian `d(pars)/d(trait)` applied on
##' top, which this function does not do.
##'
##' @return A character vector of parameter names.
##' @export
census_parameters <- function() {
  census_parameter_names()
}

## Model-specific entry point into the C++ census.
census_direct_fn <- function(strategy) {
  cls <- class(strategy)[[1]]
  switch(cls,
         FF16_Strategy = FF16_census_direct,
         TF24_Strategy = TF24_census_direct,
         TF24f_Strategy = TF24f_census_direct,
         stop("No census implemented for ", cls))
}

##' Extract a frozen stand state from a collected SCM run.
##'
##' The state is copied out as plain numeric vectors -- heights, stem number
##' densities and the two heartwood state variables, per node, per species. This
##' is what makes the direct sensitivity term "at fixed state" true by
##' construction: once the grid has left the solver, nothing a strategy or a
##' parameter assignment does can move a node's height, density or quadrature
##' weight.
##'
##' @param results Output of `run_scm(..., collect = TRUE)`.
##' @param step Time step to census. Defaults to the last.
##' @return A list with one element per species, each a data frame of `height`,
##'   `density`, `area_heartwood` and `mass_heartwood`; plus attributes `step`
##'   and `time`.
##' @export
stand_census_state <- function(results, step = NULL) {
  spp <- results$species
  if (is.null(spp) || !is.data.frame(spp)) {
    stop("`results` must be the output of run_scm(..., collect = TRUE)")
  }
  if (is.null(step)) {
    step <- max(spp$step)
  }
  at <- spp[spp$step == step, , drop = FALSE]
  if (nrow(at) == 0L) {
    stop("No state at step ", step)
  }
  labels <- unique(at$species)
  state <- lapply(labels, function(i) {
    s <- at[at$species == i, , drop = FALSE]
    data.frame(height = as.numeric(s$height),
               density = as.numeric(s$density),
               area_heartwood = as.numeric(s$area_heartwood),
               mass_heartwood = as.numeric(s$mass_heartwood))
  })
  names(state) <- paste0("species_", labels)
  attr(state, "step") <- step
  attr(state, "time") <- at$time[[1]]
  state
}

## Per-species C++ call, on a frozen state.
census_one_species <- function(state_i, strategy) {
  census_direct_fn(strategy)(strategy,
                             state_i$height,
                             state_i$density,
                             state_i$area_heartwood,
                             state_i$mass_heartwood)
}

##' Stand census on a frozen state.
##'
##' Total leaf area (`area_leaf`), stem or basal area (`area_stem`) and
##' above-ground mass (`mass_above_ground`), summed over species. Each is a
##' trapezium quadrature over the size distribution, taken on a height-sorted
##' grid: reserve-gated growth lets a younger node overtake an older one, and on
##' a crossed grid neighbouring trapezia cancel rather than accumulate.
##'
##' @param state A frozen state from `stand_census_state()`.
##' @param strategies A list of strategies, one per species, in the same order.
##' @return A named numeric vector of census totals, with a `grid_is_monotone`
##'   attribute (one logical per species) recording whether the guard was needed.
##' @export
stand_census <- function(state, strategies) {
  if (length(state) != length(strategies)) {
    stop("Need exactly one strategy per species in `state`")
  }
  parts <- Map(census_one_species, state, strategies)
  total <- Reduce(`+`, lapply(parts, "[[", "census"))
  attr(total, "grid_is_monotone") <-
    vapply(parts, "[[", logical(1), "grid_is_monotone")
  total
}

##' Direct (fixed-state) term of a stand census sensitivity.
##'
##' The derivative of each census metric with respect to each registered
##' strategy parameter, **holding the size distribution fixed**:
##'
##' \deqn{\left.\frac{\partial C}{\partial \phi}\right|_{state}
##'       = \sum_{s}\sum_k w_k n_k \frac{\partial m}{\partial \phi}}
##'
##' with \eqn{w_k} the quadrature weight, \eqn{n_k} the stem number density and
##' \eqn{m} the per-individual quantity the metric sums.
##'
##' **This is one of two terms.** The full derivative also contains the response
##' of the size distribution itself, \eqn{\sum_k \partial(w_k n_k)/\partial\phi
##' \cdot m_k}, which requires the stand's development and is not computed
##' anywhere in this package. The direct term is separable from it, is checkable
##' on its own against a finite difference on a frozen state, and is the term a
##' sweep over the size distribution can never produce.
##'
##' Columns are named `species_<label>.<parameter>`. The species prefix is not
##' decoration: without it a multi-species result repeats each parameter name
##' once per species, and character indexing then silently resolves every name
##' to species one's column.
##'
##' A registered parameter that reaches nothing is reported as an exact zero and
##' flagged `FALSE` in `support`. Where `support` is `TRUE` the accumulator ran,
##' so a zero there is a real zero. An unknown parameter is refused by name.
##'
##' @param state A frozen state from `stand_census_state()`.
##' @param strategies A list of strategies, one per species, in the same order.
##' @param parameters Optional subset of `census_parameters()` to report. Names
##'   not in the registry are an error.
##' @return A list with `census` (the metric totals), `gradient` (metrics x
##'   species-parameters) and `support` (a logical matrix of the same shape).
##' @export
stand_census_direct_term <- function(state, strategies, parameters = NULL) {
  if (length(state) != length(strategies)) {
    stop("Need exactly one strategy per species in `state`")
  }
  registered <- census_parameters()
  if (is.null(parameters)) {
    parameters <- registered
  } else {
    ## Validate against the bare registry, not against the assembled column
    ## names: the column names carry a species prefix and repeat nothing, which
    ## is exactly what makes this check able to see an unknown name at all.
    unknown <- setdiff(parameters, registered)
    if (length(unknown) > 0L) {
      stop("Unknown census parameter(s): ",
           paste(sQuote(unknown), collapse = ", "),
           ". Registered: ", paste(registered, collapse = ", "))
    }
  }

  parts <- Map(census_one_species, state, strategies)
  metrics <- census_metrics()
  labels <- names(state)
  if (is.null(labels)) {
    labels <- paste0("species_", seq_along(state))
  }

  gradient <- do.call(cbind, lapply(parts, function(x)
    x$gradient[, parameters, drop = FALSE]))
  support <- do.call(cbind, lapply(parts, function(x)
    x$support[, parameters, drop = FALSE]))
  cols <- unlist(lapply(labels, function(l) paste0(l, ".", parameters)))
  dimnames(gradient) <- list(metrics, cols)
  dimnames(support) <- list(metrics, cols)

  census <- Reduce(`+`, lapply(parts, "[[", "census"))
  list(census = census,
       gradient = gradient,
       support = support,
       grid_is_monotone = vapply(parts, "[[", logical(1), "grid_is_monotone"))
}
