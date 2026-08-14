## Three things the first pass got wrong or left open:
##  1. the dropped channel must be taken per cohort height, not only through
##     h_max -- a placement tied to cohort tops drops a term for every one;
##  2. a fixed grid's resolution is set by the canopy it is under, so it has to
##     be swept over the canopy heights a run actually passes through;
##  3. an append must be shown inert rather than assumed.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

stand <- function(n = 8, top = 18, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}
heights_of <- function(p)
  unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))

gl <- gl_nodes(220)

## ---- 1. the dropped channel, per cohort height ----------------------------
## The placement is a function of the whole height vector. Perturbing one
## cohort's height moves whichever knots that placement ties to it, with the
## profile held fixed. Summed in quadrature over cohorts, this is the whole term
## the passive positions drop.
dropped_per_cohort <- function(patch, place_of_heights, hs, z, eps = 1e-4) {
  rng <- diff(range(exact_field(patch, z, "L")$value))
  vapply(seq_along(hs), function(j) {
    hi <- hs; hi[[j]] <- hs[[j]] * (1 + eps)
    lo <- hs; lo[[j]] <- hs[[j]] * (1 - eps)
    a <- herm_eval(fit_at(patch, place_of_heights(hi), "L"), z)$value
    b <- herm_eval(fit_at(patch, place_of_heights(lo), "L"), z)$value
    max(abs(a - b)) / (2 * eps) / rng
  }, numeric(1))
}

p  <- stand()
hs <- heights_of(p)
z  <- test_grid(max(hs), hs)

places <- list(
  `uniform x h_max, 65` = function(h) place_uniform_hmax(max(h), 65),
  `breaks + fill, 65`   = function(h) place_breaks_fill(max(h), h, 65),
  `adaptive on value`   = function(h) place_adaptive(p, max(h), 1e-6, 1e-4),
  `fixed abs, d=0.25`   = function(h) place_fixed_abs(max(h), 0.25),
  `fixed abs, d=0.10`   = function(h) place_fixed_abs(max(h), 0.10)
)

cat("=== dropped position channel, per cohort height (max/range) ===\n")
cat(sprintf("%-22s %s\n", "placement",
            paste(sprintf("%9.2f", hs), collapse = "")))
for (nm in names(places)) {
  d <- dropped_per_cohort(p, places[[nm]], hs, z)
  cat(sprintf("%-22s %s   | quad-sum %.3e\n", nm,
              paste(sprintf("%9.2e", d), collapse = ""), sqrt(sum(d^2))))
}

## ---- 2. append inertness ---------------------------------------------------
## Append a knot at the top of a fixed grid and require every query at or below
## the canopy to be bit-identical. Queries never exceed h_max, so the span an
## append opens carries none of them, provided a knot already sits at or above
## the canopy.
cat("\n=== append inertness, fixed abs d=0.25 ===\n")
for (top_pad in 0:3) {
  base <- place_fixed_abs(max(hs), 0.25, pad = top_pad)
  more <- place_fixed_abs(max(hs), 0.25, pad = top_pad + 1)
  a <- herm_eval(fit_at(p, base, "L"), z)$value
  b <- herm_eval(fit_at(p, more, "L"), z)$value
  cat(sprintf("pad %d -> %d knots, top knot %6.2f (h_max %5.2f): max |change| %.3e  bit-identical %s\n",
              top_pad, length(base), max(base), max(hs),
              max(abs(a - b)), identical(a, b)))
}

## ---- 3. the canopy-height sweep -------------------------------------------
## A fixed grid's count and relative resolution both move with the canopy, so
## the question it has to answer is whether it holds across the heights a run
## passes through -- not at one of them.
cat("\n=== consumer error (max relative crown-mean) vs canopy height ===\n")
cat(sprintf("%6s %7s %6s | %-11s %-11s | %-11s %-11s | %-11s %-11s\n",
            "h_max", "L(0)", "nCoh", "unif65", "", "fixed .25", "", "fixed .10", ""))
cat(sprintf("%6s %7s %6s | %5s %5s | %5s %5s | %5s %5s | %5s %5s\n",
            "", "", "", "knots", "crown", "knots", "crown", "knots", "crown", "", ""))
for (top in c(1.5, 3, 6, 12, 18, 25, 35)) {
  pp <- stand(top = top, floor_h = min(0.6, top / 6), ldens = -1.6)
  h  <- heights_of(pp); hm <- max(h)
  eta <- ladder_strategy_parameter(pp, 1, "eta")
  zz <- test_grid(hm, h)
  ef <- function(q) exact_field(pp, q, "L")$value
  me <- vapply(h, function(x) crown_mean(ef, x, eta, gl), numeric(1))
  row <- c()
  for (pl in list(function() place_uniform_hmax(hm, 65),
                  function() place_fixed_abs(hm, 0.25),
                  function() place_fixed_abs(hm, 0.10))) {
    kx <- pl()
    f <- fit_at(pp, kx, "L")
    ff <- function(q) herm_eval(f, q)$value
    mf <- vapply(h, function(x) crown_mean(ff, x, eta, gl), numeric(1))
    row <- c(row, length(kx), max(abs(mf - me) / me))
  }
  cat(sprintf("%6.1f %7.4f %6d | %5d %9.2e | %5d %9.2e | %5d %9.2e\n",
              hm, ef(0), length(h), row[1], row[2], row[3], row[4], row[5], row[6]))
}
