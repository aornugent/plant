##' Census metrics of a solved stand.
##'
##' The trapezium integral of \eqn{n_k \psi(\mathrm{state}_k)} over the size
##' distribution, for each metric \eqn{\psi}. The grid is the coordinate the
##' density is carried in -- the birth date, or the height -- and the inflow
##' boundary node closes it, so the interval between it and the smallest cohort
##' is in the sum.
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
##' with. The cohort-height columns carry the trapezium weights as well as the
##' integrand only where the height is the coordinate integrated over; on the
##' birth-date coordinate the weights are constants.
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

##' Traits the sweep reaches no equation for, and so refuses by name.
##'
##' The leaf is entered at a solved operating point with its derivatives supplied
##' rather than recorded, so a trait of the leaf's own reaches the tape only if a
##' derivative is supplied for it. Every leaf trait the strategy declares as
##' differentiable now has one, so this set is empty. An exact zero here would be
##' indistinguishable from a trait the model genuinely does not use, which is the
##' one failure shape this design cannot afford, so a trait that loses its row
##' belongs here rather than in the matrix.
##'
##' A parameter the strategy does not declare as differentiable is a different
##' case and is already refused by name as an unknown trait.
##'
##' @return A character vector of trait names, without the species prefix a
##'   gradient's columns carry.
##' @export
stand_gradient_unanswered <- function() {
  character(0)
}

# The parameter a gradient column names, with its species index stripped. The
# refusal above is a property of the strategy's parameter and not of which
# species carries it, so it is matched on this rather than on the column name.
trait_without_species <- function(x) sub("^[0-9]+\\.", "", x)

##' Gradient of a stand's census metrics with respect to traits.
##'
##' Doubles in and doubles out: the active scalar is created and destroyed inside
##' one call and never crosses into R.
##'
##' Two terms, because a census reads the traits as well as the state. The reverse
##' pass is seeded on the states the census reads at the end of the run and runs
##' back over the steps the adaptive pass resolved; added to it is the census's own
##' reading of the traits at that state, which is not a sensitivity of the state
##' and so no sweep produces it.
##'
##' @param scm A run \code{SCM} object for the TF24 strategy.
##' @param metrics Census metric names to differentiate; defaults to all of them.
##' @param traits Column names to differentiate with respect to; defaults to
##'   every differentiable parameter every species declares. A column is named
##'   for its species and its parameter, as \code{"1.lma"}.
##' @return A list with \code{value}, the metrics at the end of the run;
##'   \code{gradient}, a metrics-by-traits matrix; and \code{control}, the four
##'   entries the gradient was taken at.
##' @export
stand_gradient <- function(scm, metrics = NULL, traits = NULL) {
  all_metrics <- census_metric_names_tf24()
  # Every column the strategy declares stays in the matrix and is CLASSIFIED,
  # rather than being dropped. A caller comparing shapes, or indexing by
  # position, must see the same width whatever the sweep can currently answer;
  # what changes is the class attached to a column, not whether it exists.
  all_traits <- census_trait_names_tf24(scm)
  if (is.null(metrics)) {
    metrics <- all_metrics
  }
  # Asking for everything is not asking for the unanswered ones: the default
  # takes what the sweep can answer and reports the rest as a class. Naming one
  # explicitly is a different question and gets a refusal.
  asked_by_name <- !is.null(traits)
  if (is.null(traits)) {
    traits <- all_traits
  }
  unknown <- setdiff(metrics, all_metrics)
  if (length(unknown) > 0L) {
    stop("Unknown census metric: ", paste(unknown, collapse = ", "))
  }
  unanswered <- traits[trait_without_species(traits) %in%
                         stand_gradient_unanswered()]
  if (asked_by_name && length(unanswered) > 0L) {
    stop("The sweep supplies no derivative for these traits, so it refuses ",
         "them rather than returning a zero that reads as a finding: ",
         paste(unanswered, collapse = ", "))
  }
  unknown <- setdiff(traits, all_traits)
  if (length(unknown) > 0L) {
    # A bare parameter name is the shape of the defect the prefix exists to
    # prevent, so say what the columns are called rather than only refusing.
    hint <- if (any(unknown %in% trait_without_species(all_traits))) {
      paste0(". Columns carry their species index, as \"",
             all_traits[[1]], "\"")
    } else {
      ""
    }
    stop("Unknown trait: ", paste(unknown, collapse = ", "), hint)
  }

  value <- stand_census(scm)[metrics]
  gradient <- do.call(rbind, census_trait_gradient_tf24(scm))
  dimnames(gradient) <- list(all_metrics, all_traits)

  list(value = value,
       gradient = gradient[metrics, traits, drop = FALSE],
       unanswered = unanswered,
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
