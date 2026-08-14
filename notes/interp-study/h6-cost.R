## What the placement costs, measured rather than argued.
##
## Two things scale with the knot count K and they are different quantities:
##  - the field build, which is O(K x cohorts) reduction evaluations per stage;
##  - the recorded block's input width, n_cohort_reads() = 2K + layers, which is
##    materialised per cohort per stage on the differentiated path.
## The forward query cost does not: a crown integral touches a fixed number of
## quadrature points whatever K is, and the grid stays equally spaced so the
## O(1) arithmetic index survives.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

stand <- function(n = 8, top = 18, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

bench <- function(f, reps) {
  f(); # warm
  t <- system.time(for (i in seq_len(reps)) f())[["elapsed"]]
  t / reps
}

cat("=== field build, one compute_environment() call ===\n")
cat(sprintf("%-22s %6s %12s %10s\n", "policy", "knots", "sec/build", "rel"))
base <- NA
settings <- list(
  list("canopy", 0.10,  65, 2), list("canopy", 0.10, 129, 2),
  list("canopy", 0.10, 257, 2),
  list("fixed",  0.25,  65, 2), list("fixed",  0.10,  65, 2),
  list("fixed",  0.05,  65, 2))
for (s in settings) {
  do.call(interp_policy_set, s)
  p <- stand()
  p$compute_environment()
  nk <- length(ladder_field_knots_tf24(p)$height)
  tm <- bench(function() p$compute_environment(), 300)
  if (is.na(base)) base <- tm
  cat(sprintf("%-22s %6d %12.3e %9.2fx\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]),
              nk, tm, tm / base))
}

## ---- the fraction of a whole run the build is ------------------------------
cat("\n=== whole SCM run ===\n")
run_once <- function() {
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(seq(0, 1.8, by = 0.15), seq(0, 1.8, by = 0.15))
  ladder_run(p)
}
cat(sprintf("%-22s %6s %12s %10s\n", "policy", "knots", "sec/run", "rel"))
base_run <- NA
for (s in list(list("canopy", 0.10, 65, 2),
               list("fixed", 0.25, 65, 2),
               list("fixed", 0.10, 65, 2))) {
  do.call(interp_policy_set, s)
  st <- try(run_once(), silent = TRUE)
  if (inherits(st, "try-error")) { cat(sprintf("%-22s   run failed\n", s[[1]])); next }
  nk <- length(ladder_field_knots_tf24(ladder_as_patch(st))$height)
  tm <- bench(function() run_once(), 3)
  if (is.na(base_run)) base_run <- tm
  cat(sprintf("%-22s %6d %12.3e %9.2fx\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]),
              nk, tm, tm / base_run))
}
