## Blast radius on the files that read the light field most directly. Every model
## shares this one construct, so a placement change re-blesses all of them; this
## is how many assertions notice.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design",
                    quiet = TRUE)
  library(testthat)
})

mode  <- Sys.getenv("INTERP_MODE", "canopy")
delta <- as.numeric(Sys.getenv("INTERP_DELTA", "0.10"))
ceil  <- as.numeric(Sys.getenv("INTERP_CEIL", "0"))
interp_policy_set(mode, delta, 65L, 1L, ceil)
cat(sprintf("=== %s delta=%.3f ceiling=%.0f ===\n", mode, delta, ceil))

d <- "/home/a/dev/plant-dev/plant/tests/testthat"
files <- c("test-environment.R", "test-canopy-methods.R",
           "test-strategy-ff16.R", "test-strategy-ff16-reference-comparison.R",
           "test-scm.R", "test-patch.R", "test-environment-TF24.R")
for (f in files) {
  t0 <- Sys.time()
  r <- tryCatch(as.data.frame(test_file(file.path(d, f), reporter = "silent")),
                error = function(e) NULL)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(r)) { cat(sprintf("%-46s ERROR (%.0fs)\n", f, el)); next }
  cat(sprintf("%-46s pass %4d  fail %3d  (%.0fs)\n", f,
              sum(r$passed), sum(r$failed) + sum(r$error), el))
}
