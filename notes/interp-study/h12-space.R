## The second design axis: the shipped field applies Beer's law BEFORE
## interpolating, so the knots hold L = exp(-A). Holding the optical depth A
## instead moves the exponential to the query. Two things could follow -- a
## different accuracy, and a positivity that is structural rather than clamped.
##
## Compared on the consumer's own metric, which is the only one commensurate
## across the two spaces.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

gl <- gl_nodes(220)

stand <- function(n = 8, top = 18.037, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

compare <- function(patch, label) {
  hs <- unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  hmax <- max(hs)
  eta <- ladder_strategy_parameter(patch, 1, "eta")
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- vapply(hs, function(h) crown_mean(ex, h, eta, gl), numeric(1))
  z <- test_grid(hmax, hs)

  cat(sprintf("\n=== %s  (h_max %.2f, L(0) %.4g) ===\n", label, hmax, ex(0)))
  cat(sprintf("%-22s %6s %12s %12s %12s\n", "placement", "knots",
              "crown (L)", "crown (A)", "min fit L"))
  for (nm in names(list(a = 1))) NULL
  cands <- list(
    `uniform x h_max, 65` = place_uniform_hmax(hmax, 65),
    `fixed abs, d=0.25`   = place_fixed_abs(hmax, 0.25),
    `fixed abs, d=0.10`   = place_fixed_abs(hmax, 0.10))
  for (nm in names(cands)) {
    kx <- cands[[nm]]
    fL <- fit_at(patch, kx, "L")
    fA <- fit_at(patch, kx, "A")
    mL <- vapply(hs, function(h) crown_mean(function(q) herm_eval(fL, q)$value,
                                            h, eta, gl), numeric(1))
    mA <- vapply(hs, function(h) crown_mean(function(q) exp(-herm_eval(fA, q)$value),
                                            h, eta, gl), numeric(1))
    cat(sprintf("%-22s %6d %12.3e %12.3e %12.3e\n", nm, length(kx),
                max(abs(mL - me) / me), max(abs(mA - me) / me),
                min(herm_eval(fL, z)$value)))
  }
}

compare(stand(),                        "shaded, L(0)=0.055")
compare(stand(ldens = -1.0),            "deeper shade")
compare(stand(ldens = -0.39),           "deepest: L(0)~6e-05")
compare(ladder_patch_two_by_two(),      "open stand")

## Does the interpolant ever go negative, which is what the clamp at
## resource_spline.h's get_value_at_height guards? A field held in A cannot,
## because exp is positive whatever the interpolant does.
cat("\n=== how close the L-space interpolant comes to zero ===\n")
for (ld in c(-1.6, -1.0, -0.39, 0.5, 1.5)) {
  p <- stand(ldens = ld)
  hs <- unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  z <- test_grid(max(hs), hs)
  f <- fit_at(p, place_fixed_abs(max(hs), 0.10), "L")
  v <- herm_eval(f, z)$value
  exv <- exact_field(p, z, "L")$value
  cat(sprintf("L(0) %10.3e   min interpolated %12.4e   min exact %12.4e   undershoot %s\n",
              exv[[1]], min(v), min(exv), min(v) < 0))
}
