# The light field's knot spacing is its only discretisation parameter, so
# refining it and watching the model output is the convergence test. Run this to
# reproduce the table in NEWS.md.
#
#   Rscript notes/light-field-convergence.R
#
# Every environment fixes its own spacing, and ResourceSpline takes it on
# construction, so a study varies it by handing the environment a different one
# before the run.

library(plant)

reference_stand <- function(spacing) {
  p0 <- scm_base_parameters("FF16")
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                       birth_rate = list(20))
  env <- Environment("FF16")
  if (!is.null(spacing)) {
    env$light_availability <- ResourceSpline(spacing)
  }
  scm <- run_scm(p1, env, Control())
  list(offspring = scm$offspring_production,
       knots = nrow(scm$patch$environment$light_availability$state))
}

spacings <- c(0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)
res <- lapply(spacings, reference_stand)
converged <- res[[length(res)]]$offspring

cat(sprintf("%-9s %-7s %-14s %s\n", "spacing", "knots", "offspring", "rel. to converged"))
for (i in seq_along(spacings)) {
  cat(sprintf("%-9.5f %-7d %-14.9f %.2e\n", spacings[i], res[[i]]$knots,
              res[[i]]$offspring,
              abs(res[[i]]$offspring - converged) / converged))
}

# The fitness landscape, which is what the field is read for. The optimum must
# not move between the shipped spacing and a refinement of it.
lma <- exp(seq(log(0.03), log(0.6), length.out = 13))
landscape <- function(spacing) {
  vapply(lma, function(l) {
    p0 <- scm_base_parameters("FF16")
    p0$max_patch_lifetime <- 40
    p1 <- add_strategies(p0, trait_matrix(l, "lma"), hyperpar = FF16_hyperpar,
                         birth_rate = list(1.0))
    env <- Environment("FF16")
    env$light_availability <- ResourceSpline(spacing)
    run_scm(p1, env, Control())$net_reproduction_ratios
  }, numeric(1))
}
shipped <- landscape(0.05)
refined <- landscape(0.0125)
cat(sprintf("\noptimum lma: shipped %.5g, refined %.5g\n",
            lma[which.max(shipped)], lma[which.max(refined)]))
cat(sprintf("worst relative difference over the landscape: %.2e\n",
            max(abs(shipped - refined) / refined)))
