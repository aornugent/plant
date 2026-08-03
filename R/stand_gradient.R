##' Census metrics of a solved stand.
##'
##' The trapezium integral of \eqn{n_k \psi(\mathrm{state}_k)} over the cohort
##' heights, for each metric \eqn{\psi}. The quadrature grid starts at the inflow
##' boundary node, so the interval between it and the smallest cohort is in the
##' sum.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A named numeric vector, one entry per census metric.
##' @export
stand_census <- function(scm) {
  stats::setNames(census_tf24(scm), census_metric_names_tf24())
}

##' Sensitivity of the census metrics to the ODE state at the current time.
##'
##' One row per metric and one column per ODE state entry, in the order
##' \code{scm$patch$ode_state} writes. This is what the reverse pass is seeded
##' with, and the cohort-height columns carry the trapezium weights as well as
##' the integrand.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A numeric matrix, metrics by ODE state entries.
##' @export
stand_census_state_adjoint <- function(scm) {
  rows <- census_state_adjoint_tf24(scm)
  out <- do.call(rbind, rows)
  rownames(out) <- census_metric_names_tf24()
  out
}

##' The four \code{Control} entries a census gradient depends on.
##'
##' Each of these moves the trajectory the gradient is taken along, so two
##' gradients taken at different values of any of them are gradients of
##' different functions. \code{stand_gradient} records them and
##' \code{stand_gradient_compare} refuses a pair that disagrees.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @return A named numeric vector of length four.
##' @export
gradient_control <- function(scm) {
  stats::setNames(gradient_control_tf24(scm),
                  c("GSS_tol_abs", "ci_abs_tol", "node_gradient_eps",
                    "schedule_eps"))
}

##' Gradient of a stand's census metrics with respect to traits.
##'
##' Doubles in and doubles out: the active scalar is created and destroyed inside
##' one call and never crosses into R. The census functional is seeded on the
##' states it reads at the end of the run, and the reverse pass runs back over
##' the steps the adaptive pass resolved.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @param metrics Census metric names to differentiate; defaults to all of them.
##' @param traits Trait names to differentiate with respect to; defaults to every
##'   differentiable parameter the strategy declares.
##' @return A list with \code{value}, the metrics at the end of the run;
##'   \code{gradient}, a metrics-by-traits matrix; and \code{control}, the four
##'   entries the gradient was taken at.
##' @export
stand_gradient <- function(scm, metrics = NULL, traits = NULL) {
  all_metrics <- census_metric_names_tf24()
  all_traits <- census_trait_names_tf24(scm)
  if (is.null(metrics)) {
    metrics <- all_metrics
  }
  if (is.null(traits)) {
    traits <- all_traits
  }
  unknown <- setdiff(metrics, all_metrics)
  if (length(unknown) > 0L) {
    stop("Unknown census metric: ", paste(unknown, collapse = ", "))
  }
  unknown <- setdiff(traits, all_traits)
  if (length(unknown) > 0L) {
    stop("Unknown trait: ", paste(unknown, collapse = ", "))
  }

  value <- stand_census(scm)[metrics]
  ## SPIKE (p3/trait-mask): `traits` now reaches C++. Only the requested trait
  ## columns are computed, and only they come back -- an unrequested column is
  ## absent, not zero.
  gradient <- do.call(rbind, census_trait_gradient_tf24(scm, traits, traits))
  dimnames(gradient) <- list(all_metrics, traits)
  if (any(!is.finite(gradient))) {
    stop("stand_gradient: a requested trait column is not finite")
  }

  list(value = value,
       gradient = gradient[metrics, , drop = FALSE],
       control = gradient_control(scm))
}

##' Compare two stand gradients.
##'
##' Two gradients are comparable only if they were taken at the same
##' \code{Control}: each of the four entries changes the trajectory and so
##' changes the function being differentiated.
##'
##' @param a,b Results of \code{stand_gradient}.
##' @return The element-wise difference \code{a$gradient - b$gradient}, over the
##'   metrics and traits both carry.
##' @export
stand_gradient_compare <- function(a, b) {
  differing <- names(a$control)[!identical_doubles(a$control, b$control)]
  if (length(differing) > 0L) {
    stop("These gradients were taken at different Control values, so they are ",
         "gradients of different functions: ",
         paste(differing, collapse = ", "))
  }
  metrics <- intersect(rownames(a$gradient), rownames(b$gradient))
  traits <- intersect(colnames(a$gradient), colnames(b$gradient))
  a$gradient[metrics, traits, drop = FALSE] -
    b$gradient[metrics, traits, drop = FALSE]
}

# Element-wise exact equality, NA-safe, for two named numeric vectors of the
# same names.
identical_doubles <- function(x, y) {
  vapply(names(x), function(n) identical(unname(x[[n]]), unname(y[[n]])),
         logical(1))
}
