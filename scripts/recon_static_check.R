# Static per-RK-stage env reconstruction check (#472 scope B, coupled-replay R0 de-risk).
#
# Question: does reconstructing competition from the HARVESTED per-RK-stage stand
# (heights h_i + per-node competition effects ce_i) at the env spline's OWN frozen
# knot x-positions, then Beer's law light = exp(-comp), reproduce the SCM's stored
# per-stage env knot VALUES to ~3e-14? The trapezium is the same arithmetic as
# Species::compute_competition -- EXCEPT the harvest (r_compute_competition_effect_by_nodes
# / r_heights) iterates `nodes` only and OMITS the boundary new_node, whose tail term
# Species::compute_competition DOES include. This isolates that gap before building
# the coupled stepper (the stepper's R0 faithfulness is bounded by it).
suppressMessages({library(devtools); load_all(".", compile = FALSE, quiet = TRUE)})

p <- scm_base_parameters("FF16")
p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                    birth_rate = list(20))
p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
scm <- run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE),
               refine_schedule = FALSE)

patch <- scm$patch
eh  <- patch$environment_history
shs <- patch$stand_height_stage_history
scs <- patch$stand_competition_stage_history
nnh <- patch$stand_newnode_height_stage_history       # boundary node height [step][stage]
nnc <- patch$stand_newnode_competition_stage_history  # boundary node effect [step][stage]
eta <- scm$parameters$strategies[[1]]$pars$eta

stopifnot(length(shs) == length(eh), length(scs) == length(eh))

# Trapezium competition at z over a stand (descending heights h, effects ce),
# matching Species::compute_competition INCLUDING the boundary new_node (hb, ceb)
# as the lowest node (its tail segment evaluates to the SCM's conditional tail term:
# above the boundary top g_b = 0, so appending it unconditionally is value-identical).
recon_comp <- function(z, h, ce, hb, ceb) {
  h <- c(h, hb); ce <- c(ce, ceb)
  ord <- order(h, decreasing = TRUE)
  h <- h[ord]; ce <- ce[ord]
  Q <- function(u) ifelse(u < 1, (1 - u^eta)^2, 0)
  g <- ce * Q(z / h)
  if (length(h) < 2) return(0)
  0.5 * sum((h[-length(h)] - h[-1]) * (g[-length(g)] + g[-1]))
}

# Walk every (step, stage), compare reconstructed knot light to the SCM's stored
# knot y-values. Report the worst absolute + relative deviation, and where it lands.
worst_abs <- 0; worst_rel <- 0; worst_loc <- ""
n_knots_total <- 0
per_stage_max <- numeric(0)
for (n in seq_along(eh)) {
  for (s in 1:6) {
    e  <- eh[[n]][[s]]
    h  <- shs[[n]][[s]]
    ce <- scs[[n]][[s]]
    if (length(h) == 0) next
    xk <- e$light_availability$spline$x
    yk <- e$light_availability$spline$y        # = SCM's exp(-compute_competition(xk))
    recon <- vapply(xk, function(z) exp(-recon_comp(z, h, ce, nnh[[n]][[s]], nnc[[n]][[s]])), numeric(1))
    ae <- abs(recon - yk)
    re <- ae / pmax(abs(yk), 1e-300)
    per_stage_max <- c(per_stage_max, max(ae))
    n_knots_total <- n_knots_total + length(xk)
    if (max(ae) > worst_abs) {
      worst_abs <- max(ae)
      j <- which.max(ae)
      worst_rel <- re[j]
      worst_loc <- sprintf("step %d stage %d knot %d (z=%.4f, yk=%.6e, recon=%.6e, nnodes=%d)",
                           n, s, j, xk[j], yk[j], recon[j], length(h))
    }
  }
}

cat(sprintf("\nStatic reconstruction check over %d steps x 6 stages (%d knots):\n",
            length(eh), n_knots_total))
cat(sprintf("  worst |recon - SCM| (abs) = %.3e\n", worst_abs))
cat(sprintf("  rel there                = %.3e\n", worst_rel))
cat(sprintf("  at: %s\n", worst_loc))
cat(sprintf("  median per-stage max-abs = %.3e\n", median(per_stage_max)))
cat(sprintf("  3e-14 gate: %s\n", if (worst_abs < 1e-12) "PASS" else "FAIL (boundary node / other gap)"))

# Where does the gap live? Recon should match best HIGH in the canopy (the boundary
# node, a seedling at ~h0, only shades z < h0). Tabulate worst-abs by knot height band.
cat("\nGap vs knot height (is it concentrated near the ground / boundary node?):\n")
bands <- list(ground = c(0, 1), mid = c(1, 5), high = c(5, Inf))
acc <- setNames(numeric(3), names(bands))
for (n in seq_along(eh)) for (s in 1:6) {
  e <- eh[[n]][[s]]; h <- shs[[n]][[s]]; ce <- scs[[n]][[s]]
  if (length(h) == 0) next
  xk <- e$light_availability$spline$x; yk <- e$light_availability$spline$y
  recon <- vapply(xk, function(z) exp(-recon_comp(z, h, ce, nnh[[n]][[s]], nnc[[n]][[s]])), numeric(1))
  ae <- abs(recon - yk)
  for (b in names(bands)) {
    sel <- xk >= bands[[b]][1] & xk < bands[[b]][2]
    if (any(sel)) acc[b] <- max(acc[b], max(ae[sel]))
  }
}
print(acc)
