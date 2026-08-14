## Can a grid track the canopy's SCALE without any position being a function of
## the state?
##
## Take one fixed lattice, z = j * dmin, and use every 2^m-th point of it, with
## the level m chosen from the canopy. Every position is then a constant of the
## run -- a member of one fixed set -- while the SPACING follows the canopy, so
## the resolution is relative like canopy-following knots and the count is
## roughly constant.
##
## The cost is that a level change is a discrete event that is NOT inert: it
## coarsens spans a query reads. That is the thing to measure, against the
## dropped channel a canopy-following grid carries at every stage.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

gl <- gl_nodes(220)

stand <- function(top, floor_h, ldens = -1.6, n = 8) {
  hs <- seq(top, floor_h, length.out = n)
  ladder_patch(species = "fast", heights = list(hs),
               log_densities = list(ldens + seq(0, 0.7, length.out = n)))
}

## Level chosen so the count is at most `target`; positions stay on the lattice.
place_dyadic <- function(hmax, dmin = 0.0125, target = 200) {
  m <- 0
  while (hmax / (dmin * 2^m) + 2 > target) m <- m + 1
  d <- dmin * 2^m
  seq(0, by = d, length.out = ceiling(hmax / d) + 2)
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

cat("=== consumer error vs canopy height, at a comparable knot budget ===\n")
cat(sprintf("%7s | %-18s | %-18s | %-18s\n",
            "h_max", "canopy x h_max 65", "fixed d=0.05", "dyadic (<=200)"))
cat(sprintf("%7s | %8s %9s | %8s %9s | %8s %9s\n",
            "", "knots", "err", "knots", "err", "knots", "err"))
for (top in c(0.7, 1.5, 3, 6, 12, 18, 35)) {
  p <- stand(top, min(0.4, top / 5))
  hm <- max(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1)))))
  a <- place_uniform_hmax(hm, 65)
  b <- place_fixed_abs(hm, 0.05)
  d <- place_dyadic(hm)
  cat(sprintf("%7.1f | %8d %9.2e | %8d %9.2e | %8d %9.2e\n",
              hm, length(a), consumer(p, a), length(b), consumer(p, b),
              length(d), consumer(p, d)))
}

## The position channel. For the dyadic lattice it must be exactly zero away
## from a level change, because every position is a constant.
cat("\n=== dropped position channel (max/range), per cohort quad-sum ===\n")
p <- stand(12, 0.4)
hs <- unlist(lapply(p$species, function(s)
  vapply(s$nodes, function(n) n$height, numeric(1))))
z <- test_grid(max(hs), hs)
rng <- diff(range(exact_field(p, z, "L")$value))
chan <- function(place) sqrt(sum(vapply(seq_along(hs), function(j) {
  hi <- hs; hi[[j]] <- hs[[j]] * 1.0001
  lo <- hs; lo[[j]] <- hs[[j]] * 0.9999
  max(abs(herm_eval(fit_at(p, place(hi), "L"), z)$value -
          herm_eval(fit_at(p, place(lo), "L"), z)$value)) / 2e-4 / rng
}, numeric(1))^2))
cat(sprintf("  canopy x h_max, 65 : %.3e\n", chan(function(h) place_uniform_hmax(max(h), 65))))
cat(sprintf("  fixed d=0.05       : %.3e\n", chan(function(h) place_fixed_abs(max(h), 0.05))))
cat(sprintf("  dyadic (<=200)     : %.3e\n", chan(function(h) place_dyadic(max(h)))))

## Where a level change sits, and how big the jump in the field is when it fires.
cat("\n=== level changes over a run, and what one costs ===\n")
hs_seq <- exp(seq(log(0.4), log(35), length.out = 400))
lev <- vapply(hs_seq, function(h) length(place_dyadic(h)), numeric(1))
sw <- hs_seq[which(diff(lev) < 0) + 1]
cat(sprintf("  canopy heights at which the level coarsens: %s\n",
            paste(sprintf("%.2f", sw), collapse = ", ")))
for (h in sw) {
  p2 <- stand(h * 1.02, min(0.4, h / 5))
  hm <- max(unlist(lapply(p2$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1)))))
  before <- place_dyadic(hm * 0.99); after <- place_dyadic(hm * 1.01)
  if (length(before) == length(after)) next
  zz <- test_grid(hm, hm)
  jump <- max(abs(herm_eval(fit_at(p2, before, "L"), zz)$value -
                  herm_eval(fit_at(p2, after, "L"), zz)$value)) /
          diff(range(exact_field(p2, zz, "L")$value))
  cat(sprintf("  at h_max %.2f: %d -> %d knots, field jump %.3e of range\n",
              hm, length(before), length(after), jump))
}
