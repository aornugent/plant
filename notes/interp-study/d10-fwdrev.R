## Forward tangent against the reverse sweep, at every state a trajectory
## reached, for whichever build PLANT_BUILD names.
##
## The rung-3 check samples six states. One of them disagrees at 4e-11 on the
## lattice and at 8e-15 on the branch it was cut from, which is either a property
## of the placement or one unlucky state. Sampling every state separates those:
## a placement effect moves the whole distribution, an unlucky state moves one
## point of it.
BUILD <- Sys.getenv("PLANT_BUILD")
stopifnot(nzchar(BUILD))
suppressMessages({ library(odelia); pkgload::load_all(BUILD, quiet = TRUE) })
source(file.path(BUILD, "tests/testthat/helper-gradient-ladder.R"))

stand <- ladder_stand_trajectory(lifetime = 2)
trajectory <- stand$store_trajectory()
patch <- ladder_as_patch(stand)
names_state <- ladder_rate_names(patch)
n <- patch$ode_size

at <- seq(2, length(trajectory))
res <- data.frame(k = integer(0), worst = numeric(0), cell = character(0),
                  knots = integer(0))
for (k in at) {
  patch$set_ode_state(trajectory[[k]]$state, trajectory[[k]]$time)
  invisible(patch$ode_rates)
  fwd <- ladder_rhs_state_jacobian_forward_tf24(patch)
  rev <- t(vapply(seq_len(n), function(i)
    ladder_rhs_adjoint_tf24(patch, replace(numeric(n), i, 1))$state, numeric(n)))
  rs <- pmax(apply(abs(fwd), 1, max), apply(abs(rev), 1, max), .Machine$double.xmin)
  d <- abs(fwd - rev) / rs
  w <- which(d == max(d), arr.ind = TRUE)[1, ]
  nk <- tryCatch(nrow(patch$environment$light_availability$state), error = function(e) NA)
  res <- rbind(res, data.frame(k = k, worst = max(d),
    cell = sprintf("%s by %s", names_state[w[[1]]], names_state[w[[2]]]),
    knots = if (is.null(nk)) NA_integer_ else nk))
}

cat(sprintf("\n=== %s ===\n", basename(BUILD)))
cat(sprintf("states %d   knots %s\n", nrow(res),
            paste(range(res$knots, na.rm = TRUE), collapse = "-")))
q <- quantile(res$worst, c(0.5, 0.9, 0.99, 1))
cat(sprintf("forward vs reverse, over every state:\n"))
cat(sprintf("  median %9.2e   p90 %9.2e   p99 %9.2e   max %9.2e\n",
            q[[1]], q[[2]], q[[3]], q[[4]]))
cat(sprintf("  states over 1e-13: %d of %d\n", sum(res$worst > 1e-13), nrow(res)))
bad <- res[order(-res$worst), ][seq_len(min(8, nrow(res))), ]
cat("  worst states:\n")
for (i in seq_len(nrow(bad)))
  cat(sprintf("    state %4d  %9.2e  %-46s knots %s\n",
              bad$k[i], bad$worst[i], bad$cell[i], bad$knots[i]))
cat("  cells appearing among the worst 20:\n")
tb <- sort(table(res[order(-res$worst), ][seq_len(min(20, nrow(res))), "cell"]), decreasing = TRUE)
for (nm in names(tb)) cat(sprintf("    %-46s %d\n", nm, tb[[nm]]))
