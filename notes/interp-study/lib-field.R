## Shared harness for the light-field interpolant study.
##
## Everything here referees against the model's own reduction,
## patch$compute_competition_and_slope, which is the function the field
## interpolates. So what is measured is interpolation error and nothing else --
## no quadrature error, no leaf, no solver.

suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design",
                    quiet = TRUE)
})
source("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design/tests/testthat/helper-gradient-ladder.R")

## ---- the interpolant, in R, matching odelia::hermite_interpolator ----------
## Same span polynomial and the same linear extension past either end, so a
## disagreement here is a disagreement with the shipped construct.
herm <- function(x, y, m) list(x = x, y = y, m = m)

herm_eval <- function(f, z) {
  x <- f$x; y <- f$y; m <- f$m
  n <- length(x)
  k <- findInterval(z, x, rightmost.closed = TRUE, all.inside = TRUE)
  h <- x[k + 1] - x[k]
  a <- y[k]; b <- y[k + 1]
  sa <- m[k] * h; sb <- m[k + 1] * h
  c1 <- sa
  c2 <- 3 * (b - a) - 2 * sa - sb
  c3 <- 2 * (a - b) + sa + sb
  t <- (z - x[k]) / h
  val <- a + t * (c1 + t * (c2 + t * c3))
  slp <- (c1 + t * (2 * c2 + 3 * t * c3)) / h
  lo <- z <= x[[1]]; hi <- z >= x[[n]]
  val[lo] <- y[[1]] + m[[1]] * (z[lo] - x[[1]]); slp[lo] <- m[[1]]
  val[hi] <- y[[n]] + m[[n]] * (z[hi] - x[[n]]); slp[hi] <- m[[n]]
  list(value = val, slope = slp)
}

## ---- the exact field ------------------------------------------------------
## Beer's law is applied before interpolation in the shipped code: the knots hold
## L = exp(-A) and dL/dz = -A' exp(-A). `space = "A"` holds the optical depth
## instead, which is the second design axis.
exact_field <- function(patch, z, space = c("L", "A")) {
  space <- match.arg(space)
  as <- vapply(z, function(zz) patch$compute_competition_and_slope(zz), numeric(2))
  A <- as[1, ]; Ap <- as[2, ]
  if (space == "A") return(list(value = A, slope = Ap))
  E <- exp(-A)
  list(value = E, slope = -Ap * E)
}

## A field fitted at `knots`, sampled from the same exact profile.
fit_at <- function(patch, knots, space = "L") {
  ex <- exact_field(patch, knots, space)
  herm(knots, ex$value, ex$slope)
}

## ---- placements -----------------------------------------------------------
place_uniform_hmax <- function(hmax, n = 65) seq(0, hmax, length.out = n)

## Cohort tops, plus fill in the widest gaps until the count is reached.
place_breaks_fill <- function(hmax, heights, n = 65) {
  x <- sort(unique(c(0, heights[heights > 0 & heights < hmax], hmax)))
  while (length(x) < n) {
    g <- diff(x)
    i <- which.max(g)
    x <- sort(c(x, x[[i]] + g[[i]] / 2))
  }
  x
}

## develop's criterion: bisect every flagged interval, keep it flagged while the
## midpoint's predicted value misses the target by more than atol AND rtol.
place_adaptive <- function(patch, hmax, atol = 1e-6, rtol = 1e-4,
                           nbase = 17, max_depth = 16, space = "L") {
  x <- seq(0, hmax, length.out = nbase)
  y <- exact_field(patch, x, space)$value
  flag <- c(FALSE, rep(TRUE, nbase - 1))     # interval i sits left of x[i]
  dx <- hmax / (nbase - 1)
  for (d in seq_len(max_depth)) {
    dx <- dx / 2
    if (!any(flag)) break
    xm <- x[which(flag)] - dx
    fitc <- herm(x, y, exact_field(patch, x, space)$slope)
    pred <- herm_eval(fitc, xm)$value
    ym <- exact_field(patch, xm, space)$value
    ok <- abs(ym - pred) < atol | abs(1 - pred / ym) < rtol
    newflag <- logical(0); newx <- numeric(0); newy <- numeric(0)
    j <- 1L
    for (i in seq_along(x)) {
      if (flag[[i]]) {
        newx <- c(newx, xm[[j]]); newy <- c(newy, ym[[j]])
        newflag <- c(newflag, !ok[[j]])
        newx <- c(newx, x[[i]]); newy <- c(newy, y[[i]])
        newflag <- c(newflag, !ok[[j]])
        j <- j + 1L
      } else {
        newx <- c(newx, x[[i]]); newy <- c(newy, y[[i]])
        newflag <- c(newflag, FALSE)
      }
    }
    x <- newx; y <- newy; flag <- newflag
  }
  x
}

## A grid of constants: uniform in absolute height, extended upward only.
## `pad` knots are kept at or above the canopy so an append never moves a span
## any query reads.
place_fixed_abs <- function(hmax, delta = 0.25, pad = 2) {
  n <- ceiling(hmax / delta) + pad
  seq(0, by = delta, length.out = n + 1)
}

## A grid of constants whose spacing is proportional to height: z_j = z0 r^j,
## with the ground added as one flat span. Scale-invariant, so the resolution a
## cohort gets near its own top does not depend on the canopy's height.
place_fixed_geom <- function(hmax, r = 1.10, z0 = 0.02, pad = 2) {
  J <- ceiling(log(hmax / z0) / log(r)) + pad
  c(0, z0 * r^(0:J))
}

## ---- metrics --------------------------------------------------------------
## A dense test grid: uniform, plus clustering at each cohort top, where the
## curvature breaks are.
test_grid <- function(hmax, heights, n = 8000) {
  z <- seq(0, hmax, length.out = n)
  near <- unlist(lapply(heights, function(h)
    h + c(-1, 1) %o% (h * c(1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.3))))
  z <- sort(unique(pmax(0, pmin(hmax, c(z, near, heights)))))
  z
}

err_stats <- function(fitted, exact) {
  e <- fitted - exact
  c(max = max(abs(e)), rms = sqrt(mean(e^2)))
}

## Crown-mean light for a cohort of height h.
##   M(h) = int_0^h L(z) q(z,h) dz,  q = 2 eta (1-u^eta) u^(eta-1) / h
## Substituting t = u^eta gives M = int_0^1 L(h t^(1/eta)) 2 (1-t) dt, so the
## weight is trivial and the sampling is concentrated at the crown top, which is
## where the field's structure is.
gl_nodes <- function(n) {
  ## Gauss-Legendre on [0,1] via the Golub-Welsch eigenproblem.
  i <- 1:(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  x <- rev(e$values); w <- 2 * rev(e$vectors[1, ])^2
  list(x = (x + 1) / 2, w = w / 2)
}

crown_mean <- function(eval_fn, h, eta, gl) {
  z <- h * gl$x^(1 / eta)
  sum(gl$w * 2 * (1 - gl$x) * eval_fn(z))
}

## ---- the dropped position channel ------------------------------------------
## Report 03 section 3.3's falsifier: hold the profile fixed, move only the grid.
## The difference is the term the passive positions drop. Normalised as a
## response to a RELATIVE change in the canopy top, as a fraction of the field's
## range, so it is dimensionless and independent of the step.
dropped_channel <- function(patch, place_fn, hmax, z, space = "L", eps = 1e-4) {
  fhi <- fit_at(patch, place_fn(hmax * (1 + eps)), space)
  flo <- fit_at(patch, place_fn(hmax * (1 - eps)), space)
  d <- (herm_eval(fhi, z)$value - herm_eval(flo, z)$value) / (2 * eps)
  rng <- diff(range(exact_field(patch, z, space)$value))
  c(max = max(abs(d)) / rng, rms = sqrt(mean(d^2)) / rng)
}
