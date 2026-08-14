## Does the reduction have to cost knots x cohorts?
##
## The field is built by evaluating, at every knot, a sum over every cohort:
##
##   A(z) = sum_j c_j (1 - (z/h_j)^eta)^2   over the cohorts with h_j >= z
##
## so the build is O(K N) kernel evaluations, and that is why a finer grid costs
## what it costs. But the kernel expands,
##
##   (1 - (z/h)^eta)^2 = 1 - 2 (z/h)^eta + (z/h)^(2 eta)
##
## and every z-dependence is then z^eta or z^(2 eta), with the cohort dependence
## in c_j, c_j h_j^-eta and c_j h_j^-2eta. Three suffix sums over cohorts sorted
## by height give every knot at once, so the build is O(N + K).
##
## Written that way it is also numerically hopeless: h^-eta at eta = 12 spans
## 1e19 over a real stand's heights, so the suffix sums annihilate their own
## small terms. Carrying the sums ALREADY SCALED by the current height removes
## that -- every term stays O(c_j) -- and the sweep downward costs one power per
## knot rather than one per knot per cohort.
##
## This checks the algebra, the conditioning and the speed, in that order,
## because the first two decide whether the third is worth having.

eta <- 12

## The reduction as the model does it: a sum over cohorts at each knot.
naive <- function(z, c, h, eta) {
  A <- numeric(length(z)); Ap <- numeric(length(z))
  for (i in seq_along(z)) {
    k <- h >= z[i]
    if (!any(k)) next
    u <- (z[i] / h[k])^eta
    A[i] <- sum(c[k] * (1 - u)^2)
    ## q carries a 1/z that is 0/0 at the crown base; the limit is zero for
    ## every eta above 1, which is the same limit the model takes.
    Ap[i] <- if (z[i] > 0) sum(c[k] * 2 * (1 - u) * (-eta * u / z[i])) else 0
  }
  list(A = A, Ap = Ap)
}

## The same reduction as one downward sweep. T1 and T2 are the scaled sums
##   T1(z) = sum_{h_j >= z} c_j (z/h_j)^eta,  T2 likewise at 2 eta,
## every term of which is at most c_j because the sum runs over h_j >= z. Moving
## down from z to z' rescales what is already accumulated and adds the cohorts
## that have just come into range, each at its own O(1) ratio.
sweep_down <- function(z, c, h, eta) {
  ord <- order(z, decreasing = TRUE)
  zz <- z[ord]
  oh <- order(h, decreasing = TRUE)
  hs <- h[oh]; cs <- c[oh]
  A <- numeric(length(zz)); Ap <- numeric(length(zz))
  S0 <- 0; T1 <- 0; T2 <- 0
  prev <- NA_real_
  j <- 1L
  for (i in seq_along(zz)) {
    zi <- zz[i]
    if (!is.na(prev) && prev > 0) {
      rr <- (zi / prev)^eta
      T1 <- T1 * rr; T2 <- T2 * rr * rr
    }
    while (j <= length(hs) && hs[j] >= zi) {
      u <- (zi / hs[j])^eta
      S0 <- S0 + cs[j]; T1 <- T1 + cs[j] * u; T2 <- T2 + cs[j] * u * u
      j <- j + 1L
    }
    A[i] <- S0 - 2 * T1 + T2
    Ap[i] <- if (zi > 0) (-2 * eta * T1 + 2 * eta * T2) / zi else 0
    prev <- zi
  }
  o <- integer(length(zz)); o[ord] <- seq_along(zz)
  list(A = A[o], Ap = Ap[o])
}

## What a cohort consumes is L = exp(-A), so a RELATIVE error in the light is an
## ABSOLUTE error in A. Reporting a relative error in A instead would report the
## canopy top, where A goes to zero and no consumer cares, as the worst place.
## The scale to beat is the rounding the naive sum itself carries, eps * A(0).
cat("=== agreement, over stands spanning the real height range ===\n")
cat(sprintf("%6s %6s %12s %12s %12s %12s\n",
            "N", "K", "A(0)", "abs err A", "eps*A(0)", "abs err A'"))
set.seed(1)
for (cfg in list(c(8, 200), c(50, 400), c(108, 700), c(300, 1000))) {
  N <- cfg[[1]]; K <- cfg[[2]]
  ## heights piled near the canopy top with a tail to the seed, as a real stand
  h <- sort(c(0.4, 17.3 * (1 - rbeta(N - 1, 0.6, 6))))
  cc <- runif(N, 1e-4, 0.2)
  z <- seq(0, max(h) * 1.02, length.out = K)
  a <- naive(z, cc, h, eta); b <- sweep_down(z, cc, h, eta)
  cat(sprintf("%6d %6d %12.4f %12.3e %12.3e %12.3e\n", N, K, a$A[[1]],
              max(abs(a$A - b$A)), .Machine$double.eps * a$A[[1]],
              max(abs(a$Ap - b$Ap))))
}

cat("\n=== what the naive suffix sums do, for contrast ===\n")
## The unscaled form, to show the conditioning claim is real rather than assumed.
naive_separable <- function(z, c, h, eta) {
  vapply(z, function(zi) {
    k <- h >= zi
    if (!any(k)) return(0)
    sum(c[k]) - 2 * zi^eta * sum(c[k] * h[k]^-eta) + zi^(2 * eta) * sum(c[k] * h[k]^(-2 * eta))
  }, numeric(1))
}
h <- sort(c(0.4, 17.3 * (1 - rbeta(107, 0.6, 6)))); cc <- runif(108, 1e-4, 0.2)
z <- seq(0.01, max(h), length.out = 400)
a <- naive(z, cc, h, eta)$A
cat(sprintf("  unscaled suffix sums : max rel err %.3e\n",
            max(abs(naive_separable(z, cc, h, eta) - a) / pmax(abs(a), 1e-300))))
cat(sprintf("  scaled downward sweep: max rel err %.3e\n",
            max(abs(sweep_down(z, cc, h, eta)$A - a) / pmax(abs(a), 1e-300))))

cat("\n=== cost, kernel evaluations ===\n")
cat(sprintf("%6s %6s %14s %12s %10s\n", "N", "K", "naive N*K", "sweep N+K", "ratio"))
for (cfg in list(c(108, 349), c(108, 695), c(300, 700))) {
  N <- cfg[[1]]; K <- cfg[[2]]
  cat(sprintf("%6d %6d %14d %12d %9.0fx\n", N, K, N * K, N + K, (N * K) / (N + K)))
}

cat("\n=== measured time, R, same work both ways ===\n")
N <- 108; K <- 700
h <- sort(c(0.4, 17.3 * (1 - rbeta(N - 1, 0.6, 6)))); cc <- runif(N, 1e-4, 0.2)
z <- seq(0, max(h), length.out = K)
t1 <- system.time(for (i in 1:20) naive(z, cc, h, eta))[["elapsed"]]
t2 <- system.time(for (i in 1:20) sweep_down(z, cc, h, eta))[["elapsed"]]
cat(sprintf("  naive %.3f s   sweep %.3f s   %.1fx\n", t1, t2, t1 / t2))
