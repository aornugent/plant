## The canopy grid's accuracy advantage and its severed channel have one cause:
## its top knot sits exactly on the tallest cohort's curvature break. Test the
## hybrid -- a grid of constants with ONE knot snapped to h_max -- against both.
##
## The field above h_max is exactly (1, 0), and both are constants, so the
## snapped knot's DATA carries no dependence on h_max at all; only its position
## does. That is one moving knot instead of all of them.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(220)

stand <- function(n = 8, top = 18.037, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

## fixed grid, with the knots strictly between the last one below h_max and
## h_max removed, then h_max inserted: keeps the spacing regular and puts a knot
## on the break.
place_snap <- function(hmax, delta = 0.10, pad = 1) {
  x <- place_fixed_abs(hmax, delta, pad)
  keep <- x < hmax - 0.25 * delta
  sort(unique(c(x[keep], hmax, hmax + delta * seq_len(pad))))
}

metrics <- function(patch, place_of_heights, label) {
  hs <- unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  hmax <- max(hs)
  eta <- ladder_strategy_parameter(patch, 1, "eta")
  z <- test_grid(hmax, hs)
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- vapply(hs, function(h) crown_mean(ex, h, eta, gl), numeric(1))

  kx <- place_of_heights(hs)
  f <- fit_at(patch, kx, "L")
  mf <- vapply(hs, function(h)
    crown_mean(function(q) herm_eval(f, q)$value, h, eta, gl), numeric(1))
  crown <- max(abs(mf - me) / me)
  ## the tallest cohort's own crown mean, which is the one FF16's offspring
  ## production is most sensitive to
  crown_top <- abs(mf[[1]] - me[[1]]) / me[[1]]

  rng <- diff(range(ex(z)))
  d <- vapply(seq_along(hs), function(j) {
    hi <- hs; hi[[j]] <- hs[[j]] * (1 + 1e-4)
    lo <- hs; lo[[j]] <- hs[[j]] * (1 - 1e-4)
    a <- herm_eval(fit_at(patch, place_of_heights(hi), "L"), z)$value
    b <- herm_eval(fit_at(patch, place_of_heights(lo), "L"), z)$value
    max(abs(a - b)) / 2e-4 / rng
  }, numeric(1))

  cat(sprintf("%-28s %6d %11.3e %11.3e %11.3e %11.3e\n", label, length(kx),
              crown, crown_top, max(d), sqrt(sum(d^2))))
  invisible(NULL)
}

for (lab_p in list(list("shaded, L(0)=0.055", stand()),
                   list("open 2x2", ladder_patch_two_by_two()))) {
  cat(sprintf("\n=== %s ===\n", lab_p[[1]]))
  cat(sprintf("%-28s %6s %11s %11s %11s %11s\n", "placement", "knots",
              "crown max", "crown top", "drop max", "drop quad"))
  p <- lab_p[[2]]
  metrics(p, function(h) place_uniform_hmax(max(h), 65),  "canopy uniform, 65")
  metrics(p, function(h) place_uniform_hmax(max(h), 129), "canopy uniform, 129")
  metrics(p, function(h) place_fixed_abs(max(h), 0.10),   "fixed abs, d=0.10")
  metrics(p, function(h) place_snap(max(h), 0.10),        "fixed + h_max snapped, 0.10")
  metrics(p, function(h) place_fixed_abs(max(h), 0.05),   "fixed abs, d=0.05")
  metrics(p, function(h) place_snap(max(h), 0.05),        "fixed + h_max snapped, 0.05")
  metrics(p, function(h) place_breaks_fill(max(h), h, 65),"breaks + fill, 65")
}
