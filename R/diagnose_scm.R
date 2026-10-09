##' A run's node error, estimated from runs at every other and every fourth of
##' its introductions, each with its invaders walked on its recording.
##'
##' Each coarser run keeps every other introduction of the one above it, both
##' ends included, so an odd count halves the spacing exactly and an even one
##' leaves a single spacing last. Each chooses its own steps at \code{ctrl}.
##' Where the error falls at the square law, as it does on the birth-date
##' coordinate, the converged answer is near
##' \code{value + (value - every_other) / 3}: \code{error} is that correction,
##' and \code{ratio}, \code{(every_other - every_fourth) /
##' (value - every_other)}, is near 4 where the estimate holds. The two coarser
##' runs, with their walks and sweeps, cost about three quarters of the run
##' again.
##'
##' @title Diagnose a run
##' @param p Parameters of one species, carrying no recorded ODE schedule.
##' @param env,ctrl,events As for \code{\link{run_scm}}. Each run's
##'   introductions replace those among \code{events}.
##' @param invaders A named list of \code{Parameters}, each walked at each run's
##'   introductions.
##' @param gradient Should each run and walk give its elasticities too, by a
##'   sweep? TF24 only.
##' @return A list of
##'   \item{quantities}{one row per run and quantity, with \code{value},
##'     \code{every_other}, \code{every_fourth}, \code{error} and
##'     \code{ratio}: \code{"log_offspring_production"}, and with
##'     \code{gradient} for each \code{parameter} its \code{"elasticity"},
##'     d log(offspring production) / d log(theta), or where theta is 0 its
##'     \code{"derivative"}, d log(offspring production) / d theta.}
##'   \item{distance}{one row per invader and parameter it moves, from the
##'     stand whose recording it walks: on \code{scale} \code{"log_ratio"},
##'     log(theta' / theta), or where theta is 0 \code{"difference"},
##'     theta' - theta.}
##'   \item{failures}{one row per run, walk or sweep that threw or was
##'     refused, with its run's introductions and the message. Its quantities
##'     are \code{NA}.}
##'   \item{introductions}{the number of introductions of each run.}
##' @export
diagnose_scm <- function(p, env = NULL, ctrl = control(), invaders = list(),
                         events = NULL, gradient = TRUE) {
  if (length(p$strategies) != 1 || length(p$ode_times) > 0) {
    stop("diagnose_scm() takes one species with no recorded ODE schedule")
  }
  if (length(invaders) > 0 &&
      (is.null(names(invaders)) || any(names(invaders) %in% c("", "stand")))) {
    stop("diagnose_scm() names each invader, none of them \"stand\"")
  }
  every_other <- function(t) t[unique(c(seq(1, length(t), 2), length(t)))]
  times <- list(p$node_schedule_times[[1]])
  times[[2]] <- every_other(times[[1]])
  times[[3]] <- every_other(times[[2]])
  if (length(times[[3]]) == length(times[[2]])) {
    stop("diagnose_scm() needs four introductions or more")
  }

  failures <- data.frame(run = character(), introductions = integer(),
                         stage = character(), message = character())
  fail <- function(run, n, stage, message) {
    failures[nrow(failures) + 1, ] <<- list(run, n, stage, message)
  }
  measure <- function(scm, q, run, n) {
    offspring <- sum(scm$offspring_production)
    out <- data.frame(quantity = "log_offspring_production",
                      parameter = NA_character_, value = log(offspring))
    if (!gradient) return(out)
    g <- tryCatch(stand_gradient(scm, metrics = "offspring_production"),
                  error = function(e) conditionMessage(e))
    refusal <- if (is.character(g)) g else g$refusal[["offspring_production"]]
    if (!is.null(refusal)) {
      fail(run, n, "sweep", paste(unlist(refusal), collapse = " "))
      return(out)
    }
    grad <- g$gradient["offspring_production", ]
    pars <- q$strategies[[1]]$pars
    names(grad) <- trait_without_species(names(grad))
    theta <- vapply(names(grad), function(n) pars[[n]], 0)
    rbind(out, data.frame(
      quantity = ifelse(theta == 0, "derivative", "elasticity"),
      parameter = names(grad),
      value = unname(ifelse(theta == 0, 1, theta) * grad / offspring)))
  }
  ## The run at each set of introductions, then each invader walked on it.
  measured <- lapply(times, function(t) {
    n <- length(t)
    pt <- p
    pt$node_schedule_times <- list(t)
    ev <- events
    if (!is.null(ev)) {
      other <- ev$type != "node_introduction"
      ev <- events(lapply(unclass(ev), `[`, other), node_introductions(pt))
    }
    scm <- tryCatch(run_scm(pt, env, ctrl, events = ev,
                            record_trajectory = gradient),
                    error = function(e) conditionMessage(e))
    if (is.character(scm)) {
      fail("stand", n, "run", scm)
      return(list())
    }
    out <- list(stand = measure(scm, pt, "stand", n))
    for (k in names(invaders)) {
      q <- invaders[[k]]
      q$node_schedule_times <- list(t)
      walked <- tryCatch({ scm$run_mutant(q); TRUE },
                         error = function(e) conditionMessage(e))
      if (isTRUE(walked)) {
        out[[k]] <- measure(scm, q, k, n)
      } else {
        fail(k, n, "walk", walked)
      }
    }
    out
  })

  key <- function(m) paste(m$quantity, m$parameter)
  quantities <- do.call(rbind, lapply(c("stand", names(invaders)), function(k) {
    v <- lapply(measured, `[[`, k)
    rows <- unique(do.call(rbind, lapply(v, `[`, c("quantity", "parameter"))))
    if (is.null(rows)) return(NULL)
    at <- function(i) {
      if (is.null(v[[i]])) return(NA_real_)
      v[[i]]$value[match(key(rows), key(v[[i]]))]
    }
    data.frame(run = k, rows, value = at(1), every_other = at(2),
               every_fourth = at(3), error = (at(1) - at(2)) / 3,
               ratio = (at(2) - at(3)) / (at(1) - at(2)), row.names = NULL)
  }))
  theta <- unlist(p$strategies[[1]]$pars)
  distance <- do.call(rbind, lapply(names(invaders), function(k) {
    moved <- unlist(invaders[[k]]$strategies[[1]]$pars)[names(theta)]
    keep <- moved != theta
    data.frame(run = k, parameter = names(theta)[keep],
               scale = ifelse(theta == 0, "difference", "log_ratio")[keep],
               distance = unname(ifelse(theta == 0, moved - theta,
                                        log(moved / theta))[keep]),
               row.names = NULL)
  }))
  list(quantities = quantities, distance = distance, failures = failures,
       introductions = lengths(times))
}
