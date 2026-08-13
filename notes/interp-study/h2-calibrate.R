## Find (a) a realistically shaded stand for the field study and (b) the darkest
## stand whose leaves still classify interior, which is what the tangent needs.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant", quiet = TRUE)
})
source("/home/a/dev/plant-dev/plant/tests/testthat/helper-gradient-ladder.R")

shaded_patch <- function(n = 8, top = 18, floor_h = 0.6, ldens = -0.39,
                         relative_reserve = 0.12) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld),
               relative_reserve = relative_reserve)
}

probe <- function(ldens, n = 8, top = 18, floor_h = 0.6) {
  p <- try(shaded_patch(n = n, top = top, floor_h = floor_h, ldens = ldens),
           silent = TRUE)
  if (inherits(p, "try-error")) return(list(ok = FALSE, msg = "patch build"))
  A0 <- p$compute_competition_and_slope(0.0)[[1]]
  jac <- try(ladder_rhs_state_jacobian_forward_tf24(p), silent = TRUE)
  ok <- !inherits(jac, "try-error")
  msg <- if (ok) "interior" else "refused"
  list(ok = ok, A0 = A0, L0 = exp(-A0), msg = msg, patch = p)
}

cat("=== shade sweep: n=8, top=18, floor=0.6 ===\n")
for (ld in c(-4.0, -3.0, -2.5, -2.0, -1.6, -1.2, -0.39)) {
  r <- probe(ld)
  cat(sprintf("ldens %6.2f   A(0)=%8.4f  L(0)=%10.3e   tangent: %s\n",
              ld, r$A0, r$L0, r$msg))
}

cat("\n=== shade sweep: floor raised to 3 m (no deep-shade seedlings) ===\n")
for (ld in c(-3.0, -2.5, -2.0, -1.6, -1.2, -0.39)) {
  r <- probe(ld, floor_h = 3)
  cat(sprintf("ldens %6.2f   A(0)=%8.4f  L(0)=%10.3e   tangent: %s\n",
              ld, r$A0, r$L0, r$msg))
}

cat("\n=== shade sweep: 4 cohorts, top 12, floor 2 ===\n")
for (ld in c(-3.0, -2.0, -1.5, -1.0, -0.39, 0.3)) {
  r <- probe(ld, n = 4, top = 12, floor_h = 2)
  cat(sprintf("ldens %6.2f   A(0)=%8.4f  L(0)=%10.3e   tangent: %s\n",
              ld, r$A0, r$L0, r$msg))
}
