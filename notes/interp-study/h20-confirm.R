## Final confirmation of the two claims the design rests on, from a clean process.
source("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design/notes/interp-study/lib-field.R")

stand <- function(top = 18.037, ldens = -2.2, n = 8, floor_h = 0.6) {
  hs <- seq(top, floor_h, length.out = n)
  ladder_patch(species = "fast", heights = list(hs),
               log_densities = list(ldens + seq(0, 0.7, length.out = n)))
}

## 1. the position channel, per cohort
p <- stand(); hs <- unlist(lapply(p$species, function(s)
  vapply(s$nodes, function(n) n$height, numeric(1))))
z <- test_grid(max(hs), hs)
rng <- diff(range(exact_field(p, z, "L")$value))
chan <- function(place) vapply(seq_along(hs), function(j) {
  hi <- hs; hi[[j]] <- hs[[j]] * 1.0001
  lo <- hs; lo[[j]] <- hs[[j]] * 0.9999
  max(abs(herm_eval(fit_at(p, place(hi), "L"), z)$value -
          herm_eval(fit_at(p, place(lo), "L"), z)$value)) / 2e-4 / rng
}, numeric(1))
c_can <- chan(function(h) place_uniform_hmax(max(h), 65))
c_fix <- chan(function(h) place_fixed_abs(max(h), 0.10))
cat("CLAIM 1  dropped position channel, per cohort (max/range)\n")
cat(sprintf("  canopy x h_max, 65 : quad-sum %.3e   nonzero on %d of %d cohorts\n",
            sqrt(sum(c_can^2)), sum(c_can > 0), length(hs)))
cat(sprintf("  fixed abs, d=0.10  : quad-sum %.3e   nonzero on %d of %d cohorts\n",
            sqrt(sum(c_fix^2)), sum(c_fix > 0), length(hs)))
cat(sprintf("  -> all exactly zero: %s\n\n", all(c_fix == 0)))

## 2. end to end, canopy off a knot
resid <- function(args) {
  do.call(interp_policy_set, args)
  q <- stand()
  J <- ladder_rhs_state_jacobian_forward_tf24(q)
  fd <- ladder_rhs_state_difference(q, rel = 1e-6)
  nm <- ladder_rate_names(q); hc <- which(nm == "height")
  r <- vapply(hc, function(j)
    max(abs(J[, j] - fd[, j])) / max(abs(fd[, j])), numeric(1))
  c(tallest = r[[1]], median_other = median(r[-1]),
    knots = length(ladder_field_knots_tf24(q)$height))
}
a <- resid(list("canopy", 0.10, 65, 1, 0))
b <- resid(list("fixed", 0.10, 65, 1, 0))
cat("CLAIM 2  tangent vs difference of the rates, tallest height column\n")
cat(sprintf("  canopy x h_max, 65 (%d knots): tallest %.3e  median other %.3e  ratio %.0fx\n",
            a[["knots"]], a[["tallest"]], a[["median_other"]],
            a[["tallest"]] / a[["median_other"]]))
cat(sprintf("  fixed abs, d=0.10  (%d knots): tallest %.3e  median other %.3e  ratio %.1fx\n",
            b[["knots"]], b[["tallest"]], b[["median_other"]],
            b[["tallest"]] / b[["median_other"]]))
