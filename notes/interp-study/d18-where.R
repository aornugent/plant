## Where does the descent disagree with the per-height reduction?
##
## It agrees to rounding on the birth-date coordinate and not on the height one,
## and the gap grows with the knot count. That is a pattern, so locate it rather
## than reason about it: print the disagreement knot by knot against where the
## cohorts actually are.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid",
                    quiet = TRUE)
})

p0 <- scm_base_parameters("FF16")
p0$max_patch_lifetime <- 10
p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                     hyperpar = FF16_hyperpar, birth_rate = list(20))
scm <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
scm$collect <- FALSE
scm$run()
p <- scm$patch
sp <- p$species[[1]]
h <- sort(sp$heights, decreasing = TRUE)

cat(sprintf("coordinate is birth date: %s\n",
            Control()$node_density_in_birth_date))
cat(sprintf("%d cohorts, heights %.4f down to %.4f\n", length(h), h[[1]],
            h[[length(h)]]))

st <- p$environment$light_availability$state
z <- st[, "height"]
swept <- -log(st[, "light_availability"])                      # the field holds exp(-A)
per_height <- vapply(z, function(q) p$compute_competition_and_slope(q)[1],
                     numeric(1))
d <- swept - per_height

cat(sprintf("\nknots %d, worst |diff| %.3e\n", length(z), max(abs(d))))
o <- order(-abs(d))[1:12]
cat(sprintf("\n%8s %14s %14s %12s %10s %s\n", "z", "swept A", "per-height A",
            "diff", "cohorts>=z", "nearest cohort top"))
for (i in sort(o)) {
  above <- sum(h >= z[i])
  near <- h[which.min(abs(h - z[i]))]
  cat(sprintf("%8.3f %14.9f %14.9f %12.3e %10d %10.4f\n",
              z[i], swept[i], per_height[i], d[i], above, near))
}

## Is the gap where cohorts are, or everywhere?
cat(sprintf("\ndiff at knots above the canopy   : %.3e\n",
            max(abs(d[z > h[[1]]]))))
cat(sprintf("diff at knots below the shortest : %.3e\n",
            max(abs(d[z < h[[length(h)]]]))))
cat(sprintf("diff inside the cohort range     : %.3e\n",
            max(abs(d[z >= h[[length(h)]] & z <= h[[1]]]))))

## The first knot at which it appears, scanning down from the canopy: a descent
## that loses something loses it once and carries it, so the onset localises the
## step that did it.
inside <- which(z <= h[[1]])
first <- inside[which(abs(d[inside]) > 1e-12)]
if (length(first)) {
  k <- max(first)
  cat(sprintf("\nscanning down from the canopy, the gap first exceeds 1e-12 at z = %.4f\n",
              z[k]))
  cat(sprintf("  cohorts at or above that height: %d, tallest %.4f\n",
              sum(h >= z[k]), h[[1]]))
  for (j in seq(max(1, k - 3), min(length(z), k + 3))) {
    cat(sprintf("    z %8.4f  diff %12.3e  cohorts>=z %4d\n",
                z[j], d[j], sum(h >= z[j])))
  }
}
