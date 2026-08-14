## Blast radius: which of the existing tests move when the placement changes.
## A placement change re-blesses every model, because all of them read this one
## field -- so this is the practical cost of the change, not a detail.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design",
                    quiet = TRUE)
  library(testthat)
})

mode <- Sys.getenv("INTERP_MODE", "canopy")
delta <- as.numeric(Sys.getenv("INTERP_DELTA", "0.10"))
interp_policy_set(mode, delta, 65L, 1L, 0)
cat(sprintf("=== policy: %s delta=%.3f ===\n", mode, delta))

d <- "/home/a/dev/plant-dev/plant/tests/testthat"
files <- setdiff(list.files(d, "^test-.*\\.[Rr]$"),
                 c("test-mutant.R", "test-strategy-tf24.R",
                   "test-strategy-tf24f.R"))
tot_f <- 0L; tot_p <- 0L; moved <- character(0)
for (f in files) {
  r <- tryCatch(as.data.frame(test_file(file.path(d, f), reporter = "silent")),
                error = function(e) NULL)
  if (is.null(r)) { moved <- c(moved, paste0(f, " [error]")); next }
  nf <- sum(r$failed) + sum(r$error)
  np <- sum(r$passed)
  tot_f <- tot_f + nf; tot_p <- tot_p + np
  if (nf > 0) moved <- c(moved, sprintf("%s (%d)", f, nf))
}
cat(sprintf("passed %d   failed/errored %d   files with failures %d\n",
            tot_p, tot_f, length(moved)))
cat(paste(moved, collapse = "\n"), "\n")
