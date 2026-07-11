##' Reverse-mode gradients of emergent stand metrics with respect to traits.
##'
##' These differentiate a stand summary (e.g. offspring production) of a run SCM
##' with respect to a focal species' traits, using odelia's reverse-mode AD over
##' the recorded trajectory. The SCM must have been run with
##' \code{control(save_RK45_cache = TRUE)} so the competitive environment is
##' recorded.
##'
##' \code{invasion_gradient} differentiates the frozen-canopy (invasion) metrics:
##' the resident's competitive landscape is held fixed while the focal species'
##' traits are perturbed. \code{stand_gradient} is the resident workflow, in
##' which the canopy re-shades with the trait (self-shading feedback).
##' \code{offspring_production_gradient} is the invasion convenience for the
##' offspring-production metric.
##'
##' @title Emergent stand gradients
##' @param scm A run \code{SCM} object (the resident).
##' @param metrics Character vector of emergent metric names (e.g.
##'   \code{"offspring_production"}).
##' @param traits Character vector of trait names (e.g. \code{"lma"}).
##' @param species Index (1-based) of the focal species whose traits vary.
##' @return \code{invasion_gradient}/\code{stand_gradient} return a list with
##'   \code{value} (the metric values) and \code{gradient} (an m x n matrix, rows
##'   = metrics, cols = traits). \code{offspring_production_gradient} returns the
##'   offspring-production gradient over the requested traits.
##' @rdname stand_gradient
##' @export
invasion_gradient <- function(scm, metrics, traits, species = 1L) {
  strategy <- extract_RcppR6_template_types(scm, "SCM")[[1]]
  invasion_gradient_cpp(scm, as.character(metrics), as.character(traits),
                        as.integer(species), strategy)
}

##' @rdname stand_gradient
##' @export
stand_gradient <- function(scm, metrics, traits, species = 1L) {
  strategy <- extract_RcppR6_template_types(scm, "SCM")[[1]]
  stand_gradient_cpp(scm, as.character(metrics), as.character(traits),
                     as.integer(species), strategy)
}

##' @rdname stand_gradient
##' @export
offspring_production_gradient <- function(scm, traits, species = 1L) {
  res <- invasion_gradient(scm, "offspring_production", traits, species)
  drop(res$gradient)
}

##' \code{birth_rate_gradient} differentiates the resident metrics with respect
##' to a focal species' birth rate on the coupled (self-shading) replay: the
##' canopy re-shades with the density the birth rate sets, so the derivative
##' carries the demographic feedback (and can flip the sign of biomass relative
##' to the identity \code{metric / birth_rate}). With \code{metrics =
##' "net_reproduction_ratio"} this is \code{dR0/db}, the plant-side derivative for
##' the equilibrium (\code{R0 = 1}) Newton solve.
##'
##' @rdname stand_gradient
##' @export
birth_rate_gradient <- function(scm, metrics, species = 1L) {
  strategy <- extract_RcppR6_template_types(scm, "SCM")[[1]]
  res <- birth_rate_gradient_cpp(scm, as.character(metrics), as.integer(species),
                                 strategy, 1.0)
  res[c("value", "gradient")]
}
