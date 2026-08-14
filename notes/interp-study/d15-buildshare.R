## How much of a right-hand side IS the field build?
##
## The separable reduction is worth building only if the build is a large part of
## the work, and that share has to be measured rather than fitted through
## whole-run timings.
##
## Two traps, both hit before this script settled. `ode_rates` is not a step: on
## its own it reads back stored rates. `set_ode_state` is not a step either: it
## rebuilds the environment and stops. The right-hand side is the pair, so that
## is what is timed, and the build is taken as a share of the pair.
##
## The saving a separable sweep offers is (N K) / (N + K), so it is a function of
## the cohort count as well as the knot count -- which is why this has to run on
## a real stand and not on a ladder fixture, where one cohort makes the two forms
## the same cost by construction.
BUILD <- Sys.getenv("PLANT_BUILD")
MODEL <- Sys.getenv("MODEL", "FF16")
stopifnot(nzchar(BUILD))
suppressMessages({ library(odelia); pkgload::load_all(BUILD, quiet = TRUE) })

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

hyper <- if (MODEL == "TF24") TF24_hyperpar else FF16_hyperpar
lifetimes <- if (MODEL == "TF24") c(4, 10) else c(4, 10, 40)

for (lt in lifetimes) {
  p0 <- scm_base_parameters(MODEL)
  p0$max_patch_lifetime <- lt
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = hyper,
                       birth_rate = list(20))
  scm <- SCM(MODEL, paste0(MODEL, "_Env"))(p1, Environment(MODEL), Control())
  scm$collect <- FALSE
  wall <- system.time(scm$run())[["elapsed"]]
  p <- scm$patch
  hs <- unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  n <- length(hs)
  nk <- nrow(p$environment$light_availability$state)
  y <- p$ode_state; tt <- p$ode_time

  t_build <- bench(function() p$compute_environment(), 200)
  t_state <- bench(function() p$set_ode_state(y, tt), 200)
  t_rates <- bench(function() invisible(p$ode_rates), 200)
  t_rhs <- t_state + t_rates
  share <- t_build / t_rhs
  saving <- (n * nk) / (n + nk)

  cat(sprintf("\n%s lifetime %2g yr : %3d cohorts, %3d knots, canopy %.2f m, run %.1f s\n",
              MODEL, lt, n, nk, max(hs), wall))
  cat(sprintf("  set_ode_state (rebuilds the field) %10.3e s\n", t_state))
  cat(sprintf("  ode_rates                          %10.3e s\n", t_rates))
  cat(sprintf("  the pair, one right-hand side      %10.3e s\n", t_rhs))
  cat(sprintf("  compute_environment alone          %10.3e s = %.0f%% of it\n",
              t_build, 100 * share))
  cat(sprintf("  reduction work now N*K = %d, swept N+K = %d, saving %.0fx\n",
              n * nk, n + nk, saving))
  cat(sprintf("  right-hand side if the build were swept: %.3f of now\n",
              1 - share + share / saving))
}
