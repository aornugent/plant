## Choosing delta. The binding constraint is the SMALLEST cohort: its shading
## transition occupies the top ~40% of its own height, so a grid uniform in
## absolute height must resolve 0.4 * h_seed, not 0.4 * h_max.
##
## delta must not be derived from the seed height: that quantity now carries a
## derivative, and a grid position built from it would drop a term in exactly the
## direction report 04 just recovered. So delta is an independent constant.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(220)

stand <- function(top, floor_h, ldens = -1.6, n = 8) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

consumer <- function(patch, kx) {
  hs <- unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  eta <- ladder_strategy_parameter(patch, 1, "eta")
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- vapply(hs, function(h) crown_mean(ex, h, eta, gl), numeric(1))
  f <- fit_at(patch, kx, "L")
  mf <- vapply(hs, function(h)
    crown_mean(function(q) herm_eval(f, q)$value, h, eta, gl), numeric(1))
  max(abs(mf - me) / me)
}

deltas <- c(0.50, 0.25, 0.10, 0.05, 0.025)
cat("=== consumer error vs delta and canopy height ===\n")
cat("(floor cohort at 0.4 m throughout: the seed's size is what binds)\n\n")
cat(sprintf("%7s | %-13s %-13s | %s\n", "h_max", "canopy 65", "canopy 257",
            paste(sprintf("%-15s", sprintf("fixed %.3f", deltas)), collapse = "")))
for (top in c(1.5, 6, 18, 35)) {
  p <- stand(top, 0.4)
  hm <- max(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1)))))
  out <- sprintf("%7.1f | %6.2e(%3d) %6.2e(%3d) |", hm,
                 consumer(p, place_uniform_hmax(hm, 65)), 65L,
                 consumer(p, place_uniform_hmax(hm, 257)), 257L)
  for (d in deltas) {
    kx <- place_fixed_abs(hm, d)
    out <- paste0(out, sprintf(" %6.2e(%4d)", consumer(p, kx), length(kx)))
  }
  cat(out, "\n")
}

cat("\n=== and with the floor cohort raised, so the height RANGE is narrow ===\n")
cat(sprintf("%7s %7s | %-13s | %s\n", "h_max", "floor", "canopy 65",
            paste(sprintf("%-15s", sprintf("fixed %.3f", deltas)), collapse = "")))
for (fl in c(0.4, 2, 8, 14)) {
  p <- stand(18.037, fl)
  hm <- 18.037
  out <- sprintf("%7.1f %7.1f | %6.2e(%3d) |", hm, fl,
                 consumer(p, place_uniform_hmax(hm, 65)), 65L)
  for (d in deltas) {
    kx <- place_fixed_abs(hm, d)
    out <- paste0(out, sprintf(" %6.2e(%4d)", consumer(p, kx), length(kx)))
  }
  cat(out, "\n")
}
