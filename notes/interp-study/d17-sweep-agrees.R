## The field the sweep builds against the field the reduction defines.
##
## The build now descends the grid once; the per-height reduction is untouched
## and is what the transpose, R and the tests referee against. They are the same
## function, so they must agree to the per-height sum's own rounding -- and the
## scale to beat is eps * A(0), not eps, because A is a sum of that size.
##
## Checked on both models and across stand sizes, since the sweep's admission
## order is by height and the two coordinates order the quadrature differently.
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid",
                    quiet = TRUE)
})

check <- function(model, hyper, lifetime) {
  p0 <- scm_base_parameters(model)
  p0$max_patch_lifetime <- lifetime
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = hyper,
                       birth_rate = list(20))
  scm <- SCM(model, paste0(model, "_Env"))(p1, Environment(model), Control())
  scm$collect <- FALSE
  wall <- system.time(scm$run())[["elapsed"]]
  p <- scm$patch
  # The inflow boundary node is not ODE state and advances after the last field
  # build, so a field left by the run is not the field the reduction answers for
  # now. Rebuild, then compare the two at the same state.
  p$compute_environment()

  st <- p$environment$light_availability$state
  z <- st[, 1]; value <- st[, 2]; slope <- st[, 3]
  n <- length(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(q) q$height, numeric(1)))))

  ## The same field, one height at a time.
  as <- vapply(z, function(q) p$compute_competition_and_slope(q), numeric(2))
  ref_value <- exp(-as[1, ])
  ref_slope <- -(as[2, ] * ref_value)
  A0 <- as[1, 1]

  cat(sprintf("\n%s lifetime %g: %d cohorts, %d knots, canopy %.2f m, run %.2f s\n",
              model, lifetime, n, length(z), max(z), wall))
  cat(sprintf("  A(0) = %.4f, so the per-height sum's own rounding is %.2e\n",
              A0, .Machine$double.eps * A0))
  cat(sprintf("  value: max abs diff %.3e   (field range %.4f)\n",
              max(abs(value - ref_value)), diff(range(ref_value))))
  cat(sprintf("  slope: max abs diff %.3e   (slope range %.4f)\n",
              max(abs(slope - ref_slope)), diff(range(ref_slope))))
  invisible(NULL)
}

check("FF16", FF16_hyperpar, 4)
check("FF16", FF16_hyperpar, 10)
check("FF16", FF16_hyperpar, 40)
check("TF24", TF24_hyperpar, 4)
