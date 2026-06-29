# Alignment probe for the coupled replay: does the number of harvested nodes at
# (step n, stage s) equal the number of cohorts with birth_step <= n? And does the
# harvested per-stage stand match a per-cohort independent replay at theta0? This
# tells the coupled stepper exactly which cohorts are alive at each step.
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

h <- plant:::ff16_harvest(scm, 1L, NULL)
patch <- scm$patch
shs <- patch$stand_height_stage_history
birth <- h$birth_step                 # 0-based step index per cohort

cat(sprintf("N steps = %d, n cohorts = %d\n", length(shs), length(birth)))
cat(sprintf("birth_step range: [%d, %d]\n", min(birth), max(birth)))

# For each step n (1-based R index = step n-1 0-based), count harvested nodes at
# each stage vs cohorts with birth <= (n-1).
mism <- 0
for (n in seq_along(shs)) {
  n0 <- n - 1L                         # 0-based step
  nalive <- sum(birth <= n0)
  for (s in 1:6) {
    nnode <- length(shs[[n]][[s]])
    if (nnode != nalive) {
      mism <- mism + 1
      if (mism <= 8)
        cat(sprintf("  step %d(0-based %d) stage %d: harvested %d nodes, alive cohorts %d\n",
                    n, n0, s, nnode, nalive))
    }
  }
}
cat(sprintf("stage-count mismatches: %d (of %d stage-slots)\n", mism, length(shs) * 6))
