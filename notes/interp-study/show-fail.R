suppressMessages({ library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid", quiet=TRUE)
  library(testthat) })
d <- "/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid/tests/testthat"
for (f in strsplit(Sys.getenv("FILES"), ",")[[1]]) {
  res <- test_file(file.path(d,f), reporter="silent")
  for (ctx in res) for (r in ctx$results)
    if (inherits(r, "expectation_failure") || inherits(r, "expectation_error"))
      cat(sprintf("\n--- %s :: %s ---\n%s\n", f, ctx$test, substr(conditionMessage(r), 1, 700)))
}
