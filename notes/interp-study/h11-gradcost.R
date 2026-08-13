## The gradient side of the knot count. n_cohort_reads() = 2K + layers, so the
## recorded block's input vector widens with K and is materialised per cohort per
## stage. This is the one place a fixed grid's larger K could bite.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

bench <- function(f, reps) { f(); system.time(for (i in seq_len(reps)) f())[["elapsed"]] / reps }

stand <- function(n = 8, top = 18.037, floor_h = 0.6, ldens = -2.2) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

cat("=== one cohort block: input width and sweep cost ===\n")
cat(sprintf("%-20s %6s %8s %12s %8s %12s %8s\n", "policy", "knots", "inputs",
            "sec/reverse", "rel", "sec/forward", "rel"))
b1 <- NA; b2 <- NA
for (s in list(list("canopy", 0.10, 65, 1), list("canopy", 0.10, 129, 1),
               list("fixed", 0.50, 65, 1), list("fixed", 0.25, 65, 1),
               list("fixed", 0.10, 65, 1), list("fixed", 0.05, 65, 1))) {
  do.call(interp_policy_set, s)
  p <- stand()
  nk <- length(ladder_field_knots_tf24(p)$height)
  nin <- length(ladder_block_input_names_tf24(p, 1L))
  tr <- bench(function() ladder_block_jacobian_reverse_tf24(p, 1L), 40)
  tf <- bench(function() ladder_block_jacobian_forward_tf24(p, 1L), 40)
  if (is.na(b1)) { b1 <- tr; b2 <- tf }
  cat(sprintf("%-20s %6d %8d %12.3e %7.2fx %12.3e %7.2fx\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]),
              nk, nin, tr, tr / b1, tf, tf / b2))
}

cat("\n=== a whole trajectory sweep ===\n")
run_mid <- function(lifetime = 8) {
  p <- ladder_parameters(c("fast", "slow"), lifetime = lifetime)
  p$node_schedule_times <- list(seq(0, lifetime * 0.9, length.out = 10),
                                seq(0, lifetime * 0.9, length.out = 10))
  ladder_run(p)
}
cat(sprintf("%-20s %6s %12s %8s\n", "policy", "knots", "sec/gradient", "rel"))
b3 <- NA
for (s in list(list("canopy", 0.10, 65, 1), list("fixed", 0.25, 65, 1),
               list("fixed", 0.10, 65, 1))) {
  do.call(interp_policy_set, s)
  st <- run_mid()
  nk <- length(ladder_field_knots_tf24(ladder_as_patch(st))$height)
  g <- try(stand_gradient(st), silent = TRUE)
  if (inherits(g, "try-error")) {
    cat(sprintf("%-20s %6d   refused: %s\n", s[[1]], nk,
                substr(conditionMessage(attr(g, "condition")), 1, 60)))
    next
  }
  tm <- bench(function() stand_gradient(st), 2)
  if (is.na(b3)) b3 <- tm
  cat(sprintf("%-20s %6d %12.3e %7.2fx\n",
              sprintf("%s d=%.2f n=%d", s[[1]], s[[2]], s[[3]]), nk, tm, tm / b3))
}
