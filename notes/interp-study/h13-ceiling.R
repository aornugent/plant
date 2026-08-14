## The constant-K variant: a grid spanning [0, ceiling] whatever the canopy is
## doing. Positions are constants AND the count is one too, which is what report
## 05 section 5 relies on. The question is what the knots above the canopy cost.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

stand <- function(n = 8, top = 18.037, floor_h = 0.6, ldens = -2.2) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

cat("=== the knot count over a run, append-only against a ceiling ===\n")
run_long <- function(lifetime = 20) {
  p <- ladder_parameters(c("fast", "slow"), lifetime = lifetime)
  p$node_schedule_times <- list(seq(0, lifetime * 0.9, length.out = 20),
                                seq(0, lifetime * 0.9, length.out = 20))
  ladder_run(p)
}
for (s in list(list("fixed", 0.10, 65, 1, 0), list("fixed", 0.10, 65, 1, 40))) {
  do.call(interp_policy_set, s)
  st <- run_long()
  pa <- ladder_as_patch(st)
  nk <- length(ladder_field_knots_tf24(pa)$height)
  tm <- bench(function() run_long(), 3)
  cat(sprintf("ceiling %5.1f  final knots %5d  sec/run %.3e  leaf area %.7e\n",
              s[[5]], nk, tm, ladder_census_leaf_area(st)))
}

cat("\n=== does the ceiling change any number? ===\n")
interp_policy_set("fixed", 0.10, 65, 1, 0)
p1 <- stand(); r1 <- p1$ode_rates
k1 <- length(ladder_field_knots_tf24(p1)$height)
interp_policy_set("fixed", 0.10, 65, 1, 40)
p2 <- stand(); r2 <- p2$ode_rates
k2 <- length(ladder_field_knots_tf24(p2)$height)
cat(sprintf("append-only %d knots vs ceiling-40 %d knots\n", k1, k2))
cat(sprintf("rates: max |diff| %.3e   bit-identical %s\n",
            max(abs(r1 - r2)), identical(r1, r2)))

cat("\n=== cost of the knots above the canopy ===\n")
cat(sprintf("%-24s %6s %8s %12s %8s\n", "policy", "knots", "inputs",
            "sec/build", "rel"))
b <- NA
for (s in list(list("canopy", 0.10, 65, 1, 0),
               list("fixed", 0.10, 65, 1, 0),
               list("fixed", 0.10, 65, 1, 25),
               list("fixed", 0.10, 65, 1, 40),
               list("fixed", 0.10, 65, 1, 60))) {
  do.call(interp_policy_set, s)
  p <- stand()
  p$compute_environment()
  nk <- length(ladder_field_knots_tf24(p)$height)
  nin <- length(ladder_block_input_names_tf24(p, 1L))
  tm <- bench(function() p$compute_environment(), 2000)
  if (is.na(b)) b <- tm
  cat(sprintf("%-24s %6d %8d %12.3e %7.2fx\n",
              sprintf("%s d=%.2f ceil=%.0f", s[[1]], s[[2]], s[[5]]),
              nk, nin, tm, tm / b))
}

## ---- the re-blessing scope: every model shares this field -------------------
cat("\n=== FF16 shifts too: the field is shared ===\n")
ff16_run <- function() {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 20
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 18, length.out = 20))
  run_scm(p, use_ode_times = FALSE)
}
interp_policy_set("canopy", 0.10, 65, 1, 0)
a <- try(ff16_run(), silent = TRUE)
if (inherits(a, "try-error")) {
  cat("FF16 run failed:", conditionMessage(attr(a, "condition")), "\n")
} else {
  ref <- a$offspring_production
  cat(sprintf("canopy 65      offspring production %.10e\n", ref))
  for (s in list(list("canopy", 0.10, 257, 1, 0), list("fixed", 0.25, 65, 1, 0),
                 list("fixed", 0.10, 65, 1, 0))) {
    do.call(interp_policy_set, s)
    b2 <- ff16_run()
    cat(sprintf("%-14s offspring production %.10e   shift %+.3e\n",
                sprintf("%s d=%.2f", s[[1]], s[[2]]),
                b2$offspring_production, b2$offspring_production / ref - 1))
  }
}
