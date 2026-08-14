## Where does a finer lattice's cost actually go? If it is the field build, the
## build is O(knots x cohorts) and the reduction has structure worth using. If it
## is elsewhere, a grid change cannot recover it.
Sys.setenv(PLANT_BUILD = "/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design")
source("/home/a/.claude/jobs/e02c60e6/tmp/d1-lib.R")

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

## A mature FF16 stand, held fixed, so only the field build varies.
mature <- function() ff16_scm(40, collect = FALSE)$patch
interp_policy_set("canopy", 0.10, 65, 1, 0)
p0 <- mature()
y <- p0$ode_state; t <- p0$ode_time
cat(sprintf("mature stand: %d cohorts, canopy %.2f m\n",
            length(heights_of(p0)), max(heights_of(p0))))

cat(sprintf("\n%-20s %7s %13s %9s %13s %9s\n",
            "placement", "knots", "sec/build", "rel", "sec/run(40y)", "rel"))
b1 <- NA; b2 <- NA
for (s in list(list("canopy", 0.10,  65, 1, 0), list("canopy", 0.10, 257, 1, 0),
               list("fixed", 0.100, 65, 1, 0), list("fixed", 0.050, 65, 1, 0),
               list("fixed", 0.025, 65, 1, 0))) {
  do.call(interp_policy_set, s)
  p <- mature()                      # same stand, rebuilt under this policy
  p$set_ode_state(y, t)
  nk <- nrow(p$environment$light_availability$state)
  tb <- bench(function() p$compute_environment(), 400)
  tr <- bench(function() ff16_scm(40, collect = FALSE), 2)
  if (is.na(b1)) { b1 <- tb; b2 <- tr }
  lab <- if (s[[1]] == "canopy") sprintf("canopy n=%d", s[[3]]) else sprintf("fixed d=%.3f", s[[2]])
  cat(sprintf("%-20s %7d %13.3e %8.2fx %13.3e %8.2fx\n", lab, nk, tb, tb / b1, tr, tr / b2))
}

## How much of the reduction is wasted? A cohort's crown shape is 1 to within
## 1e-8 below a fifth of its own height and 0 above it, so only knots inside
## [0.2 h, h] need the full evaluation; the rest is a running constant.
hs <- sort(heights_of(p0), decreasing = TRUE)
hs <- hs[is.finite(hs) & hs > 0]
cat(sprintf("\nband-limited reduction, this stand (%d cohorts):\n", length(hs)))
for (d in c(0.10, 0.05, 0.025)) {
  K <- ceiling(max(hs) / d) + 2
  z <- seq(0, by = d, length.out = K)
  naive <- K * length(hs)
  banded <- sum(vapply(hs, function(h) sum(z > 0.2 * h & z <= h), numeric(1)))
  cat(sprintf("  d=%.3f  %5d knots   naive %9.0f   banded %9.0f   saving %5.1fx\n",
              d, K, naive, banded, naive / banded))
}
