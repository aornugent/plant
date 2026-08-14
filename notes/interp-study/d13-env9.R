## Why do the tangent and the sweep disagree on one accumulator row?
##
## The rung-3 trajectory check reports up to 1.4e-10 between two exact paths at
## four or five states of a hundred, always in the same row, on both the lattice
## and the branch it was cut from. Two exact derivatives of one evaluation
## disagree only through the arithmetic, so the question is which quantity is
## being formed as a small difference of large ones.
##
## Report 00 fact 2 and its free row for root-mediated redistribution name the
## candidate: the accumulator's rate is the SUM over layers of a signed uptake,
## so a cohort drawing from one layer and releasing into another contributes a
## total that is a cancellation. This prints the evidence at the worst state
## rather than arguing from the shape.
BUILD <- Sys.getenv("PLANT_BUILD")
stopifnot(nzchar(BUILD))
suppressMessages({ library(odelia); pkgload::load_all(BUILD, quiet = TRUE) })
source(file.path(BUILD, "tests/testthat/helper-gradient-ladder.R"))

stand <- ladder_stand_trajectory(lifetime = 2)
trajectory <- stand$store_trajectory()
patch <- ladder_as_patch(stand)
nm <- ladder_rate_names(patch)
n <- patch$ode_size

## Locate the worst state rather than assuming which one it is.
worst <- list(k = NA, v = -1)
for (k in seq(2, length(trajectory))) {
  patch$set_ode_state(trajectory[[k]]$state, trajectory[[k]]$time)
  invisible(patch$ode_rates)
  fwd <- ladder_rhs_state_jacobian_forward_tf24(patch)
  rev <- t(vapply(seq_len(n), function(i)
    ladder_rhs_adjoint_tf24(patch, replace(numeric(n), i, 1))$state, numeric(n)))
  rs <- pmax(apply(abs(fwd), 1, max), apply(abs(rev), 1, max), .Machine$double.xmin)
  d <- max(abs(fwd - rev) / rs)
  if (d > worst$v) worst <- list(k = k, v = d, fwd = fwd, rev = rev, rs = rs)
}
k <- worst$k
cat(sprintf("worst state %d of %d, disagreement %.3e\n", k, length(trajectory), worst$v))

patch$set_ode_state(trajectory[[k]]$state, trajectory[[k]]$time)
invisible(patch$ode_rates)
fwd <- worst$fwd; rev <- worst$rev
d <- abs(fwd - rev) / worst$rs
w <- which(d == max(d), arr.ind = TRUE)[1, ]
ri <- w[[1]]; ci <- w[[2]]
cat(sprintf("cell: %s by %s   (row %d, col %d)\n", nm[ri], nm[ci], ri, ci))
cat(sprintf("  forward %.17g\n  reverse %.17g\n", fwd[ri, ci], rev[ri, ci]))
cat(sprintf("  absolute difference %.3e ; row scale %.3e ; relative to the cell %.3e\n",
            abs(fwd[ri, ci] - rev[ri, ci]), worst$rs[[ri]],
            abs(fwd[ri, ci] - rev[ri, ci]) / max(abs(fwd[ri, ci]), 1e-300)))

## Is the cell itself a cancellation? Compare it against the row it sits in.
cat(sprintf("\nrow %s: max |entry| %.3e, this entry %.3e, ratio %.3e\n",
            nm[ri], worst$rs[[ri]], abs(fwd[ri, ci]),
            abs(fwd[ri, ci]) / worst$rs[[ri]]))

## The accumulator's rate is a sum over layers of a signed per-cohort draw. If
## some layers are negative the total is a difference, and that is visible in
## the model rather than inferred.
cat("\nper-layer uptake by cohort at this state (the accumulator's summands):\n")
for (si in seq_along(patch$species)) {
  s <- patch$species[[si]]
  for (j in seq_along(s$nodes)) {
    cr <- tryCatch(vapply(seq_len(9), function(i)
      s$nodes[[j]]$consumption_rate(i - 1L), numeric(1)), error = function(e) NULL)
    if (is.null(cr)) { cat("  consumption_rate not reachable from R\n"); break }
    pos <- sum(cr[cr > 0]); neg <- sum(cr[cr < 0])
    cat(sprintf("  species %d node %d h=%7.3f : sum %+.6e  |pos| %.3e |neg| %.3e  cancel %.2e\n",
                si, j, s$nodes[[j]]$height, sum(cr), pos, abs(neg),
                abs(sum(cr)) / max(pos, 1e-300)))
  }
}

## The soil states, so the positivity guard's active set is visible too.
cat("\nsoil state and rates at this step:\n")
y <- patch$ode_state; r <- patch$ode_rates
env <- tail(seq_len(n), 9)
for (i in seq_along(env))
  cat(sprintf("  %-16s state %+.6e  rate %+.6e\n", nm[env[i]], y[env[i]], r[env[i]]))
