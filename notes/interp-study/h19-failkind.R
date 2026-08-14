suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design",
                    quiet = TRUE)
  library(testthat)
})
interp_policy_set("fixed", 0.10, 65L, 1L, 40)
d <- "/home/a/dev/plant-dev/plant/tests/testthat"
for (f in c("test-canopy-methods.R", "test-strategy-ff16.R")) {
  r <- as.data.frame(test_file(file.path(d, f), reporter = "silent"))
  bad <- r[r$failed + r$error > 0, ]
  cat(sprintf("\n=== %s: %d contexts with failures ===\n", f, nrow(bad)))
  print(utils::head(bad[, c("test", "failed", "error")], 12), row.names = FALSE)
}
