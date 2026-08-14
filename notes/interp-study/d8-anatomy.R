## Where does the field's error come from, and does the dropped position channel
## shrink with knot count?
##
## Two questions the design turns on and neither has been decomposed:
##
##   1. Report 03 section 2 says A(z) is C1 but not C2, with a curvature break at
##      every distinct cohort top, and that a value-chasing refiner clusters
##      around the breaks rather than landing on them. If the breaks dominate the
##      error, the lever is WHERE knots go, not HOW MANY. If they do not, the
##      lever is count and every placement argument is about cost alone.
##
##   2. Report 03 section 3.3 makes the passive-position treatment conditional on
##      the dropped channel shrinking with knot density, and records that this has
##      never been checked. Checking it decides whether a canopy-tied grid can be
##      rescued by refinement.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(220)

stand <- function(top, floor_h, ldens = -1.6, n = 8) {
  hs <- seq(top, floor_h, length.out = n)
  ladder_patch(species = "fast", heights = list(hs),
               log_densities = list(ldens + seq(0, 0.7, length.out = n)))
}

heights_of <- function(p) unlist(lapply(p$species, function(s)
  vapply(s$nodes, function(n) n$height, numeric(1))))

## ---- 1. error per span, split by whether the span straddles a cohort top ----
## The span is the unit because the Hermite polynomial is determined span by
## span: a break inside one degrades that span and no other.
cat("=== error by span, uniform grid, canopy 12 m, 8 cohorts ===\n")
p <- stand(12, 0.4)
hs <- sort(heights_of(p))
cat(sprintf("cohort tops: %s\n\n", paste(sprintf("%.2f", hs), collapse = " ")))

cat(sprintf("%8s %7s %9s | %-28s | %-28s\n", "spacing", "knots", "spans",
            "spans WITH a cohort top", "spans WITHOUT"))
cat(sprintf("%8s %7s %9s | %5s %10s %10s | %5s %10s %10s\n",
            "", "", "", "n", "max err", "mean err", "n", "max err", "mean err"))
for (d in c(0.4, 0.2, 0.1, 0.05, 0.025)) {
  x <- place_fixed_abs(max(hs), d)
  f <- fit_at(p, x, "L")
  ## sample each span densely and take that span's own worst error
  per <- vapply(seq_len(length(x) - 1), function(k) {
    zz <- seq(x[k], x[k + 1], length.out = 40)
    max(abs(herm_eval(f, zz)$value - exact_field(p, zz, "L")$value))
  }, numeric(1))
  has <- vapply(seq_len(length(x) - 1), function(k)
    any(hs > x[k] & hs < x[k + 1]), logical(1))
  ## only spans below the canopy carry any signal at all
  live <- x[-length(x)] < max(hs)
  a <- per[has & live]; b <- per[!has & live]
  cat(sprintf("%8.3f %7d %9d | %5d %10.2e %10.2e | %5d %10.2e %10.2e\n",
              d, length(x), sum(live),
              length(a), if (length(a)) max(a) else NA, if (length(a)) mean(a) else NA,
              length(b), if (length(b)) max(b) else NA, if (length(b)) mean(b) else NA))
}

## The convergence rate separates the two populations: a smooth span should show
## order 4, a span holding a curvature break order 3. That is the claim.
cat("\n=== observed order of convergence, by span population ===\n")
ords <- function(sel) {
  ds <- c(0.4, 0.2, 0.1, 0.05, 0.025)
  e <- vapply(ds, function(d) {
    x <- place_fixed_abs(max(hs), d)
    f <- fit_at(p, x, "L")
    per <- vapply(seq_len(length(x) - 1), function(k) {
      zz <- seq(x[k], x[k + 1], length.out = 40)
      max(abs(herm_eval(f, zz)$value - exact_field(p, zz, "L")$value))
    }, numeric(1))
    has <- vapply(seq_len(length(x) - 1), function(k)
      any(hs > x[k] & hs < x[k + 1]), logical(1))
    live <- x[-length(x)] < max(hs)
    v <- if (sel) per[has & live] else per[!has & live]
    if (length(v)) max(v) else NA_real_
  }, numeric(1))
  cat(sprintf("  %-22s %s\n", if (sel) "with a cohort top" else "without",
              paste(sprintf("%9.2e", e), collapse = "")))
  cat(sprintf("  %-22s %s\n", "  -> order",
              paste(c("         ", sprintf("%9.2f", -diff(log(e)) / log(2))), collapse = "")))
}
cat(sprintf("  %-22s %s\n", "spacing", paste(sprintf("%9.3f", c(0.4,0.2,0.1,0.05,0.025)), collapse = "")))
ords(TRUE); ords(FALSE)

## ---- 2. does the dropped position channel shrink with knot count? ----------
## Held profile, moved grid. If this falls like the field error, the passive
## treatment is a discretisation error and refinement rescues a canopy-tied grid.
## If it plateaus, it is a floor and no count rescues it.
cat("\n=== dropped position channel vs knot count (canopy-tied grid) ===\n")
z <- test_grid(max(hs), hs)
cat(sprintf("%7s %12s %12s %14s\n", "knots", "channel max", "field err", "ratio ch/err"))
for (n in c(33, 65, 129, 257, 513, 1025)) {
  ch <- dropped_channel(p, function(hm) place_uniform_hmax(hm, n), max(hs), z)
  f <- fit_at(p, place_uniform_hmax(max(hs), n), "L")
  ex <- exact_field(p, z, "L")$value
  er <- max(abs(herm_eval(f, z)$value - ex)) / diff(range(ex))
  cat(sprintf("%7d %12.3e %12.3e %14.2f\n", n, ch[["max"]], er, ch[["max"]] / er))
}

cat("\n=== the same, for a lattice of constants (must be identically zero) ===\n")
for (d in c(0.1, 0.05, 0.025)) {
  ch <- dropped_channel(p, function(hm) place_fixed_abs(hm, d), max(hs), z)
  cat(sprintf("  spacing %.3f : channel max %.3e\n", d, ch[["max"]]))
}

## ---- 3. what knots ON the breaks buy, at matched count ---------------------
## If the breaks dominate, a grid holding every cohort top should beat a uniform
## grid of the same size by the order gap, not by a constant.
cat("\n=== knots on the breaks vs uniform, matched count, consumer metric ===\n")
consumer <- function(patch, kx, hts) {
  eta <- ladder_strategy_parameter(patch, 1, "eta")
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- vapply(hts, function(h) crown_mean(ex, h, eta, gl), numeric(1))
  f <- fit_at(patch, kx, "L")
  mf <- vapply(hts, function(h)
    crown_mean(function(q) herm_eval(f, q)$value, h, eta, gl), numeric(1))
  max(abs(mf - me) / me)
}
cat(sprintf("%7s | %12s %12s %10s\n", "knots", "uniform", "on breaks", "gain"))
for (n in c(17, 33, 65, 129)) {
  u <- consumer(p, place_uniform_hmax(max(hs), n), hs)
  b <- consumer(p, place_breaks_fill(max(hs), hs, n), hs)
  cat(sprintf("%7d | %12.3e %12.3e %9.1fx\n", n, u, b, u / b))
}
