## What a finer lattice costs, at the canopy a real run reaches.
Sys.setenv(PLANT_BUILD = "/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design")
source("/home/a/.claude/jobs/e02c60e6/tmp/d1-lib.R")
source("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design/tests/testthat/helper-gradient-ladder.R")

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

cat("=== FF16, 40-year run ===\n")
cat(sprintf("%-22s %7s %11s %8s %14s\n", "placement", "knots", "sec/run", "rel", "offspring"))
base <- NA
for (s in list(list("canopy", 0.10,  65, 1, 0), list("fixed", 0.100, 65, 1, 0),
               list("fixed", 0.050, 65, 1, 0), list("fixed", 0.025, 65, 1, 0),
               list("fixed", 0.0125, 65, 1, 0))) {
  do.call(interp_policy_set, s)
  scm <- ff16_scm(40, collect = FALSE)
  nk <- nrow(scm$patch$environment$light_availability$state)
  tm <- bench(function() ff16_scm(40, collect = FALSE), 2)
  if (is.na(base)) base <- tm
  lab <- if (s[[1]] == "canopy") sprintf("canopy n=%d", s[[3]]) else sprintf("fixed d=%.4f", s[[2]])
  cat(sprintf("%-22s %7d %11.3e %7.2fx %14.7e\n", lab, nk, tm, tm / base,
              scm$offspring_production))
}

## The gradient side: one cohort's recording, which is linear in the knot count.
cat("\n=== TF24: one cohort's recording, and the block's input width ===\n")
cat(sprintf("%-22s %7s %8s %12s %8s\n", "placement", "knots", "inputs", "recording", "rel"))
b <- NA
for (s in list(list("canopy", 0.10, 65, 1, 0), list("fixed", 0.100, 65, 1, 0),
               list("fixed", 0.050, 65, 1, 0), list("fixed", 0.025, 65, 1, 0))) {
  do.call(interp_policy_set, s)
  p <- ladder_patch(species = "fast", heights = list(c(18.037, 12, 6, 2)),
                    log_densities = list(c(-2.2, -1.9, -1.6, -1.3)))
  nk <- length(ladder_field_knots_tf24(p)$height)
  nin <- length(ladder_block_input_names_tf24(p, 1L))
  a <- ladder_rhs_adjoint_tf24(p, rep(1, p$ode_size))
  if (is.na(b)) b <- a$block_recording_size
  lab <- if (s[[1]] == "canopy") sprintf("canopy n=%d", s[[3]]) else sprintf("fixed d=%.4f", s[[2]])
  cat(sprintf("%-22s %7d %8d %12.0f %7.2fx\n", lab, nk, nin,
              a$block_recording_size, a$block_recording_size / b))
}
