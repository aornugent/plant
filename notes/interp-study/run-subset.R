suppressMessages({ library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid", quiet=TRUE)
  library(testthat) })
d <- "/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid/tests/testthat"
files <- if (nzchar(Sys.getenv("FILES"))) strsplit(Sys.getenv("FILES"), ",")[[1]] else
  c("test-environment.R","test-canopy-methods.R","test-strategy-ff16.R",
    "test-strategy-ff16-reference-comparison.R","test-scm.R","test-patch.R",
    "test-environment-TF24.R")
for (f in files) {
  r <- tryCatch(as.data.frame(test_file(file.path(d,f), reporter="silent")), error=function(e) e)
  if (inherits(r,"error")) { cat(sprintf("%-46s ERROR %s\n", f, substr(conditionMessage(r),1,50))); next }
  bad <- r[r$failed + r$error > 0, ]
  cat(sprintf("%-46s pass %4d fail %3d\n", f, sum(r$passed), sum(r$failed)+sum(r$error)))
  if (nrow(bad)) for (i in seq_len(nrow(bad))) cat(sprintf("      - %s (%d)\n", bad$test[i], bad$failed[i]+bad$error[i]))
}
