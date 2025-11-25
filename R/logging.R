## Standard fields:
##
## * time    (done in logr)
## * level   (done in logr)
## * package (done here)
## * pid     (done here)
## * routine {equilibrium, schedule, ...}
## * stage   (start, converged, etc.)
plant_log_info <- function(...) {
  logr::log_info(..., package="plant", pid=Sys.getpid())
}

plant_log_debug <- function(...) {
  logr::log_print(..., package="plant", pid=Sys.getpid())
}

plant_log_eq <- function(...) {
  plant_log_info(..., routine="equilibrium")
}

plant_log_viable <- function(...) {
  plant_log_info(..., routine="viable")
}

plant_log_inviable <- function(...) {
  plant_log_info(..., routine="inviable")
}

make_plant_format_log_entry <- function(colour) {
  if (colour) {
    col_routine <- crayon::yellow
    col_time    <- crayon::silver
  } else {
    col_routine <- col_time <- identity
  }
  function(event) {
    if (is.null(event$routine)) {
      msg <- event$message
    } else {
      msg <- sprintf("%s> %s", col_routine(event$routine), event$message)
    }
    sprintf("[%s] %s",
            col_time(format(event$time, "%Y-%m-%d %H:%M:%OS3")),
            msg)
  }
}

##' Activate logging with logr
##'
##' By default plant prints little information about its progress.
##' This can be modified by enabling logging.  A formatter that is
##' different to the default \code{logr::log_open} formatter is
##' selected here; it will print additional information that plant's
##' internal logging functions record.
##'
##' "Schedule" events (splitting) are sent to the DEBUG stream,
##' everything else is sent to INFO.  All events have a "routine"
##' field added to them, which is useful if sent to a Redis server
##' (using \code{logr.redis}).
##' @title Activate logging with logr
##' @param ... Additional parameters passed to \code{logr::log_open},
##' but \emph{not} \code{file_name} which is hard coded here to
##' \code{"console"}.
##' @param .message,.warning,.error Include messages, warnings or
##' errors?  By default (and in contrast to \code{logr::log_open}
##' these are disabled here.
##' @export
plant_log_console <- function(...,
                              .message=FALSE,
                              .warning=FALSE,
                              .error=FALSE) {
  logr::log_open(file_name="console", ...,
                  .message=.message,
                  .warning=.warning,
                  .error=.error,
                  .formatter=make_plant_format_log_entry(TRUE))
}

## There would be no need for this though; we'd just pass through directly:
## plant_log_redis <- function(con, key, ...) {
##   logr.redis::log_redis(con, key, ...)
## }
