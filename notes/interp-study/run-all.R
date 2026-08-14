suppressMessages({ library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid", quiet=TRUE)
  library(testthat) })
d <- "/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid/tests/testthat"
files <- list.files(d, "^test-.*\\.[Rr]$")
tp <- 0; tf <- 0
for (f in files) {
  r <- tryCatch(as.data.frame(test_file(file.path(d,f), reporter="silent")), error=function(e) e)
  if (inherits(r,"error")) { cat(sprintf("%-46s ERROR %s\n", f, substr(conditionMessage(r),1,60))); tf <- tf+1; next }
  n <- sum(r$failed)+sum(r$error); tp <- tp + sum(r$passed); tf <- tf + n
  if (n) { cat(sprintf("%-46s pass %4d fail %3d\n", f, sum(r$passed), n))
    bad <- r[r$failed+r$error>0,]; for (i in seq_len(nrow(bad))) cat(sprintf("      - %s (%d)\n", bad$test[i], bad$failed[i]+bad$error[i])) }
}
cat(sprintf("\nTOTAL pass %d  fail %d  over %d files\n", tp, tf, length(files)))
