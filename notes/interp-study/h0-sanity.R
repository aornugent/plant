## Sanity: load the worktree build, make a stand, read the exact reduction.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant", quiet = TRUE)
})
source("/home/a/dev/plant-dev/plant/tests/testthat/helper-gradient-ladder.R")

p1 <- ladder_patch_one()
cat("=== ladder_patch_one ===\n")
cat("nodes      :", length(p1$species[[1]]$nodes), "\n")
hs <- vapply(p1$species[[1]]$nodes, function(n) n$height, numeric(1))
cat("heights    :", format(hs, digits = 6), "\n")
cat("height_max :", format(max(hs), digits = 8), "\n")
A0 <- p1$compute_competition_and_slope(0.0)
cat("A(0), A'(0):", format(A0, digits = 8), "\n")
cat("L(0)       :", format(exp(-A0[[1]]), digits = 8), "\n")

p4 <- ladder_patch_two_by_two()
cat("\n=== ladder_patch_two_by_two ===\n")
hs4 <- unlist(lapply(p4$species, function(s)
  vapply(s$nodes, function(n) n$height, numeric(1))))
cat("heights    :", format(sort(hs4, decreasing = TRUE), digits = 6), "\n")
A0 <- p4$compute_competition_and_slope(0.0)
cat("A(0)       :", format(A0[[1]], digits = 8), "  L(0):",
    format(exp(-A0[[1]]), digits = 8), "\n")

## What does the actual field grid look like?
k <- ladder_field_knots_tf24(p4)
cat("\nfield knots: n =", length(k$height), " range [",
    format(min(k$height), digits = 4), ",",
    format(max(k$height), digits = 6), "]\n")
cat("first 5 x  :", format(head(k$height, 5), digits = 6), "\n")
cat("names      :", paste(names(k), collapse = ", "), "\n")
