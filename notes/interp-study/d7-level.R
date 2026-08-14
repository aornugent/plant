## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))
gl <- gl_nodes(220)
stand <- function(top, floor_h, ldens = -1.6, n = 8) {
  hs <- seq(top, floor_h, length.out = n)
  ladder_patch(species = "fast", heights = list(hs),
               log_densities = list(ldens + seq(0, 0.7, length.out = n)))
}
place_dyadic <- function(hmax, dmin = 0.0125, target = 200) {
  m <- 0
  while (hmax / (dmin * 2^m) + 2 > target) m <- m + 1
  d <- dmin * 2^m
  list(x = seq(0, by = d, length.out = ceiling(hmax / d) + 2), d = d)
}
## Find canopy heights where the SPACING doubles.
hh <- exp(seq(log(0.5), log(40), length.out = 3000))
dd <- vapply(hh, function(h) place_dyadic(h)$d, numeric(1))
sw <- hh[which(diff(dd) > 0) + 1]
cat("spacing doubles at canopy heights:", sprintf("%.3f", sw), "\n\n")
cat(sprintf("%9s %8s %8s %12s %12s %12s\n",
            "h_max", "d before", "d after", "err before", "err after", "field jump"))
for (h in sw) {
  p <- stand(h, min(0.4, h / 5))
  hm <- max(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1)))))
  b <- place_dyadic(hm * 0.999); a <- place_dyadic(hm * 1.001)
  if (b$d == a$d) { cat(sprintf("%9.3f   (no level change at this stand's h_max)\n", hm)); next }
  z <- test_grid(hm, hm)
  ex <- exact_field(p, z, "L")$value
  vb <- herm_eval(fit_at(p, b$x, "L"), z)$value
  va <- herm_eval(fit_at(p, a$x, "L"), z)$value
  rng <- diff(range(ex))
  cat(sprintf("%9.3f %8.4f %8.4f %12.3e %12.3e %12.3e\n", hm, b$d, a$d,
              max(abs(vb - ex)) / rng, max(abs(va - ex)) / rng,
              max(abs(vb - va)) / rng))
}
