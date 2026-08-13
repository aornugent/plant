## End to end: does making the positions constants remove the severed column?
##
## The referee is a plain-double difference of the rates, which re-runs the
## field build and so carries the knot-position chain the tangent drops. The
## tallest cohort's height column is the one that chain reaches.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

stand <- function(n = 8, top = 18, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

height_cols <- function(patch) which(ladder_rate_names(patch) == "height")

col_resid <- function(patch, rel = 1e-6) {
  J  <- ladder_rhs_state_jacobian_forward_tf24(patch)
  fd <- ladder_rhs_state_difference(patch, rel = rel)
  list(J = J, fd = fd,
       resid = vapply(seq_len(ncol(J)), function(j) {
         a <- J[, j]; b <- fd[, j]
         max(abs(a - b)) / max(max(abs(b)), 1e-300)
       }, numeric(1)))
}

## ---- where the tangent is finite, and how shaded that is -------------------
cat("=== finiteness of the tangent vs shade ===\n")
interp_policy_set("canopy", 0.10, 65, 2)
for (ld in c(-3.0, -2.6, -2.2, -2.0, -1.8, -1.6)) {
  p <- stand(ldens = ld)
  minlight <- min(ladder_light_at_reads(p))
  J <- try(ladder_rhs_state_jacobian_forward_tf24(p), silent = TRUE)
  fin <- if (inherits(J, "try-error")) "refused" else
         sprintf("%d/%d NaN", sum(!is.finite(J)), length(J))
  cat(sprintf("ldens %5.2f  L(0) %8.4f  min read light %7.4f   %s\n",
              ld, exact_field(p, 0, "L")$value, minlight, fin))
}

## ---- the comparison, at the darkest stand with a finite tangent ------------
run_policy <- function(tag, args, ldens) {
  do.call(interp_policy_set, args)
  p <- stand(ldens = ldens)
  nk <- length(ladder_field_knots_tf24(p)$height)
  r <- col_resid(p)
  hc <- height_cols(p)
  cat(sprintf("\n--- %s  (%d knots, L(0)=%.4f) ---\n", tag, nk,
              exact_field(p, 0, "L")$value))
  cat(sprintf("%6s %14s %14s %12s\n", "col", "tangent", "difference", "rel resid"))
  for (i in seq_along(hc)) {
    j <- hc[[i]]
    k <- which.max(abs(r$fd[, j]))
    cat(sprintf("%6d %14.6e %14.6e %12.3e%s\n", j, r$J[k, j], r$fd[k, j],
                r$resid[[j]], if (i == 1) "   <- tallest" else ""))
  }
  oth <- setdiff(seq_len(ncol(r$J)), hc)
  cat(sprintf("worst non-height column: %.3e\n", max(r$resid[oth])))
  invisible(list(tallest = r$resid[[hc[[1]]]],
                 others = median(r$resid[hc[-1]])))
}

LD <- -2.2
a <- run_policy("canopy-following, 65 knots", list("canopy", 0.10, 65, 2), LD)
b <- run_policy("fixed absolute, d = 0.10",   list("fixed", 0.10, 65, 2), LD)
c_ <- run_policy("fixed absolute, d = 0.25",  list("fixed", 0.25, 65, 2), LD)

cat("\n=== summary: the tallest cohort's height column ===\n")
cat(sprintf("%-28s %12s %12s %10s\n", "placement", "tallest", "median other", "ratio"))
for (nm in list(list("canopy-following, 65", a), list("fixed abs d=0.10", b),
                list("fixed abs d=0.25", c_))) {
  cat(sprintf("%-28s %12.3e %12.3e %9.1fx\n", nm[[1]],
              nm[[2]]$tallest, nm[[2]]$others,
              nm[[2]]$tallest / nm[[2]]$others))
}

## ---- the NaN's structure, since it bounds what can be refereed in shade ----
cat("\n=== structure of the shaded-stand non-finite tangent ===\n")
interp_policy_set("canopy", 0.10, 65, 2)
for (ld in c(-2.2, -2.0, -1.6)) {
  p <- stand(ldens = ld)
  J <- ladder_rhs_state_jacobian_forward_tf24(p)
  nm <- ladder_rate_names(p)
  bad_col <- apply(J, 2, function(v) sum(!is.finite(v)))
  bad_row <- apply(J, 1, function(v) sum(!is.finite(v)))
  cat(sprintf("\nldens %.2f  min read light %.4f  NaN %d/%d\n", ld,
              min(ladder_light_at_reads(p)), sum(!is.finite(J)), length(J)))
  cat(sprintf("  columns wholly NaN: %d   partly: %d   clean: %d\n",
              sum(bad_col == nrow(J)), sum(bad_col > 0 & bad_col < nrow(J)),
              sum(bad_col == 0)))
  cat(sprintf("  rows wholly NaN   : %d   partly: %d   clean: %d\n",
              sum(bad_row == ncol(J)), sum(bad_row > 0 & bad_row < ncol(J)),
              sum(bad_row == 0)))
  aff <- unique(nm[bad_row > 0])
  cat(sprintf("  rate rows affected : %s\n", paste(aff, collapse = ", ")))
  cat(sprintf("  state cols affected: %s\n",
              paste(unique(nm[bad_col > 0]), collapse = ", ")))
}
