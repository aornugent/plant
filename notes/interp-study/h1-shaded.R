## Build a shaded stand and reproduce the tallest-cohort height-column defect.
##
## The channel is field-mediated, so it has no power on an open stand: L ~ 1
## everywhere means a crown integral cannot feel the grid move. Shade first,
## then measure.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant", quiet = TRUE)
})
source("/home/a/dev/plant-dev/plant/tests/testthat/helper-gradient-ladder.R")

## A one-species stand of many cohorts, dense enough to close the canopy.
## Heights descend so the reduction's ordering invariant holds; densities are
## distinct so a scatter reaching the wrong node shows.
shaded_patch <- function(n = 8, top = 18, floor_h = 0.6, ldens = 2.3,
                         relative_reserve = 0.12) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)   # distinct, non-commensurate-ish
  ladder_patch(species = "fast",
               heights = list(hs),
               log_densities = list(ld),
               relative_reserve = relative_reserve)
}

report <- function(patch, label) {
  hs <- vapply(patch$species[[1]]$nodes, function(n) n$height, numeric(1))
  A0 <- patch$compute_competition_and_slope(0.0)[[1]]
  cat(sprintf("%-22s n=%2d  h_max=%7.3f  A(0)=%7.4f  L(0)=%9.3e\n",
              label, length(hs), max(hs), A0, exp(-A0)))
  invisible(NULL)
}

cat("=== canopy closure vs density ===\n")
for (ld in c(-4.0, -3.0, -2.0, -1.6)) {
  report(shaded_patch(ldens = ld), sprintf("log_density base %.2f", ld))
}

## ---- the column-1 measurement --------------------------------------------
## Forward tangent state Jacobian against a plain-double difference of the
## rates. Every height column should agree; the tallest cohort's is the one
## carrying the severed knot-position chain.
column_residuals <- function(patch, rel = 1e-6) {
  J  <- ladder_rhs_state_jacobian_forward_tf24(patch)
  fd <- ladder_rhs_state_difference(patch, rel = rel)
  nm <- ladder_rate_names(patch)
  n  <- ncol(J)
  scale <- max(abs(fd))
  out <- data.frame(col = seq_len(n), name = nm,
                    jac = NA_real_, fd = NA_real_, resid = NA_real_)
  for (j in seq_len(n)) {
    a <- J[, j]; b <- fd[, j]
    den <- max(max(abs(b)), scale * 1e-12)
    out$jac[[j]]   <- a[which.max(abs(b))]
    out$fd[[j]]    <- b[which.max(abs(b))]
    out$resid[[j]] <- max(abs(a - b)) / den
  }
  out
}

p <- shaded_patch(ldens = -1.6)
report(p, "measurement stand")
res <- column_residuals(p)

hgt <- res[res$name == "height", ]
cat("\n=== height columns, tangent vs difference of the rates ===\n")
print(hgt, row.names = FALSE, digits = 4)

oth <- res[res$name != "height", ]
cat(sprintf("\nworst non-height column: %s  resid %.3e\n",
            oth$name[which.max(oth$resid)], max(oth$resid)))
cat(sprintf("tallest-cohort height column resid: %.4e\n", hgt$resid[[1]]))
cat(sprintf("median other height column resid  : %.4e\n",
            median(hgt$resid[-1])))
cat(sprintf("ratio                             : %.1f x\n",
            hgt$resid[[1]] / median(hgt$resid[-1])))
