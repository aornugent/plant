## The interpolant itself, against the kind develop uses.
##
## Every placement measurement so far fitted a HERMITE interpolant at every
## candidate knot set -- develop's adaptive criterion included -- so what those
## established is that develop's PLACEMENT loses to a lattice. The interpolant
## change was never measured, and a develop-facing case for the forward model
## rests on exactly that.
##
## Two schemes at the SAME knots, so nothing but the interpolant differs:
##
##   hermite       value AND slope at each knot, both exact, from one construct
##   solved-slope  values only; slopes solved so the result is C2 (splinefun)
##
## Hermite reads 2n numbers off the model where the solved cubic reads n, and the
## extra n are exact because the crown density is the exact vertical derivative of
## the cumulative profile. So the question is not whether the value survives -- it
## does for both -- but how much that second half buys per knot, and whether the
## slope is purchasable with knots in one scheme and not the other.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid",
                    quiet = TRUE)
})

## ---- the two interpolants, on one knot set ---------------------------------
herm_eval <- function(x, y, m, z) {
  n <- length(x)
  k <- findInterval(z, x, rightmost.closed = TRUE, all.inside = TRUE)
  h <- x[k + 1] - x[k]
  a <- y[k]; b <- y[k + 1]
  sa <- m[k] * h; sb <- m[k + 1] * h
  c1 <- sa
  c2 <- 3 * (b - a) - 2 * sa - sb
  c3 <- 2 * (a - b) + sa + sb
  t <- (z - x[k]) / h
  list(value = a + t * (c1 + t * (c2 + t * c3)),
       slope = (c1 + t * (2 * c2 + 3 * t * c3)) / h)
}

## Values only, slopes solved for C2 continuity -- the scheme develop's field is
## built on, and the one report 03 section 3.1 says has no local support.
solved_eval <- function(x, y, z) {
  f <- stats::splinefun(x, y, method = "fmm")
  list(value = f(z), slope = f(z, deriv = 1))
}

## ---- the field, and the crown-mean a cohort reads --------------------------
gl_nodes <- function(n) {
  i <- 1:(n - 1); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  xx <- rev(e$values); w <- 2 * rev(e$vectors[1, ])^2
  list(x = (xx + 1) / 2, w = w / 2)
}
GL <- gl_nodes(240)
crown_mean <- function(Lfun, h, eta) sum(GL$w * 2 * (1 - GL$x) * Lfun(h * GL$x^(1 / eta)))

run <- function(lifetime) {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- lifetime
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                       hyperpar = FF16_hyperpar, birth_rate = list(20))
  scm <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
  scm$collect <- FALSE; scm$run()
  list(patch = scm$patch, eta = p1$strategies[[1]]$pars[["eta"]])
}

exact <- function(p, z) {
  as <- vapply(z, function(q) p$compute_competition_and_slope(q), numeric(2))
  E <- exp(-as[1, ])
  list(value = E, slope = -(as[2, ] * E))
}

for (lt in c(4, 10, 40)) {
  r <- run(lt); p <- r$patch; eta <- r$eta
  hs <- sort(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(q) q$height, numeric(1)))))
  hs <- hs[is.finite(hs) & hs > 0]
  hmax <- max(hs)
  ## A dense probe grid, clustered at the cohort tops where the curvature breaks
  ## are, because that is where the two schemes are claimed to part company.
  zz <- sort(unique(pmax(0, pmin(hmax, c(seq(0, hmax, length.out = 4000),
    unlist(lapply(hs, function(h) h * (1 + c(-0.03, -0.01, -0.003, 0, 0.003, 0.01, 0.03)))))))))
  ex <- exact(p, zz)
  rng <- diff(range(ex$value))
  srng <- diff(range(ex$slope))

  cat(sprintf("\n=== FF16 lifetime %g: %d cohorts, canopy %.2f m ===\n",
              lt, length(hs), hmax))
  cat(sprintf("%7s %7s | %11s %11s %7s | %11s %11s %7s | %11s %11s\n",
              "spacing", "knots", "value herm", "value solv", "gain",
              "slope herm", "slope solv", "gain", "crown herm", "crown solv"))
  for (d in c(0.4, 0.2, 0.1, 0.05)) {
    kx <- seq(0, by = d, length.out = ceiling(hmax / d) + 3)
    ek <- exact(p, kx)
    h1 <- herm_eval(kx, ek$value, ek$slope, zz)
    s1 <- solved_eval(kx, ek$value, zz)
    ve <- c(max(abs(h1$value - ex$value)), max(abs(s1$value - ex$value))) / rng
    se <- c(max(abs(h1$slope - ex$slope)), max(abs(s1$slope - ex$slope))) / srng
    ## What a cohort actually consumes: the crown-mean of the light.
    cm <- vapply(c("h", "s"), function(which) {
      f <- if (which == "h")
        function(q) herm_eval(kx, ek$value, ek$slope, q)$value
      else function(q) solved_eval(kx, ek$value, q)$value
      me <- vapply(hs, function(h) crown_mean(function(q) exact(p, q)$value, h, eta), numeric(1))
      mf <- vapply(hs, function(h) crown_mean(f, h, eta), numeric(1))
      max(abs(mf - me) / me)
    }, numeric(1))
    cat(sprintf("%7.3f %7d | %11.2e %11.2e %6.1fx | %11.2e %11.2e %6.1fx | %11.2e %11.2e\n",
                d, length(kx), ve[1], ve[2], ve[2] / ve[1],
                se[1], se[2], se[2] / se[1], cm[1], cm[2]))
  }
}
