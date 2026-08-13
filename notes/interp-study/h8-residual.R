## Why the fixed grid's tallest column stops at 3e-05 instead of the 1.5e-06
## floor. Two candidates, and they predict opposite things:
##   (a) the perturbation straddles ceil(h_max/delta), so the knot COUNT changes
##       across the difference. Then choosing h_max off a knot removes it.
##   (b) residual interpolation error varying with the state. Then it scales with
##       delta and cannot be removed by moving h_max.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

stand <- function(n = 8, top = 18, floor_h = 0.6, ldens = -2.2) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

resid_of <- function(patch, rel = 1e-6) {
  J  <- ladder_rhs_state_jacobian_forward_tf24(patch)
  fd <- ladder_rhs_state_difference(patch, rel = rel)
  nm <- ladder_rate_names(patch)
  r <- vapply(seq_len(ncol(J)), function(j)
    max(abs(J[, j] - fd[, j])) / max(max(abs(fd[, j])), 1e-300), numeric(1))
  list(resid = r, names = nm, hc = which(nm == "height"))
}

cat("=== does h_max sitting on a knot explain it? ===\n")
cat(sprintf("%-10s %7s %8s %8s %12s %12s\n",
            "delta", "h_max", "h/delta", "straddle", "tallest", "median other"))
for (d in c(0.10, 0.25)) {
  for (top in c(18.0, 18.037, 18.5, 18.617)) {
    interp_policy_set("fixed", d, 65, 2)
    p <- stand(top = top)
    q <- top / d
    ## a relative step of 1e-6 on h_max moves it by this much
    step <- max(abs(top), 1) * 1e-6
    straddle <- ceiling((top - step) / d) != ceiling((top + step) / d)
    r <- resid_of(p)
    cat(sprintf("%-10.2f %7.3f %8.2f %8s %12.3e %12.3e\n", d, top, q,
                straddle, r$resid[[r$hc[[1]]]], median(r$resid[r$hc[-1]])))
  }
}

cat("\n=== what the worst non-height column is ===\n")
interp_policy_set("fixed", 0.10, 65, 2)
p <- stand(top = 18.037)
r <- resid_of(p)
oth <- setdiff(seq_len(length(r$resid)), r$hc)
ord <- oth[order(r$resid[oth], decreasing = TRUE)][1:6]
for (j in ord)
  cat(sprintf("col %3d  %-40s %12.3e\n", j, r$names[[j]], r$resid[[j]]))
