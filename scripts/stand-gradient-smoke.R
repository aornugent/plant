# The whole-run census trait gradient at a short patch lifetime, and the re-run
# central difference it is compared against.
#
#   Rscript scripts/stand-gradient-smoke.R <plant worktree> [lifetime]
#
# The lifetime is short so the run is cheap. This is a smoke test of the
# adjoint, not V4: V4 is the production lifetime against a re-run difference,
# and that reference is recorded as not converged in its step.
#
# The sweep crosses node introductions: it runs one segment per state width and
# contracts each newcomer's adjoint rows through the inflow boundary condition
# between segments. The driver prints the recorded widths and the segment count
# beside the gradient, since a run of one width would exercise none of it.
args <- commandArgs(TRUE)
tree <- args[1]
if (is.na(tree)) stop("give the plant worktree as the first argument")
lifetime <- if (length(args) > 1) as.numeric(args[2]) else 5

suppressMessages({ library(odelia); pkgload::load_all(tree, quiet = TRUE, export_all = TRUE) })
Sys.setenv(TESTTHAT_PARALLEL = "false")

# One trait per channel: lma reaches leaf area through the allometry, k_I
# reaches every metric through the canopy alone, and hmat through the height at
# maturity. A gradient row of exact zeros for one of these is the failure mode
# the census design is most exposed to, so the traits are named rather than
# swept.
smoke_traits <- c("lma", "k_I", "hmat")
smoke_fd_eps <- 1e-5
smoke_intro_times <- function(lifetime) seq(0, lifetime * 0.5, length.out = 8)

smoke_scm <- function(lifetime, trait = NULL, value = NULL) {
  ctrl <- Control()
  s <- TF24_Strategy()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  if (!is.null(trait)) {
    pars <- s$pars
    pars[[trait]] <- value
    s$pars <- pars
  }
  p <- scm_base_parameters("TF24")
  p$strategies <- list(s)
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list(smoke_intro_times(lifetime))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), ctrl)
  scm$run()
  scm
}

peak_mb <- function() {
  kb <- grep("^VmHWM", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", kb)) / 1024
}

scm <- smoke_scm(lifetime)
cat(sprintf("lifetime %g   nodes %d   ode_size %d   time %.6f\n",
            lifetime, length(scm$patch$species[[1]]$nodes),
            scm$patch$ode_size, scm$patch$time))
cat("\ncensus value\n")
print(stand_census(scm))
# R0 is not a census metric: a census is the trapezium integral of
# n_k psi(state_k), and R0 is a density-weighted sum of a node's own fecundity
# state, which is a different reduction. Printed as a value, not differentiated.
cat("\nR0 (value only; not in the census codomain)\n")
print(scm$net_reproduction_ratios)
cat("\ncontrol\n")
print(gradient_control(scm))

widths <- vapply(scm$store_trajectory(), function(z) length(z$state), numeric(1))
runs <- rle(widths)
cat(sprintf("\nrecorded state widths: %s   (steps: %s)   system width %d\n",
            paste(runs$values, collapse = " "), paste(runs$lengths, collapse = " "),
            scm$patch$ode_size))
cat(sprintf("width changes the sweep crosses: %d\n", length(runs$values) - 1L))

before_mb <- peak_mb()
g <- tryCatch(stand_gradient(scm, traits = smoke_traits), error = identity)
after_mb <- peak_mb()
if (inherits(g, "error")) {
  cat("\nstand_gradient refused: ", conditionMessage(g), "\n", sep = "")
} else {
  cat("\ngradient (metrics x traits)\n")
  print(g$gradient)
  # A row of exact zeros is what a severed accumulator, a snapshot read and an
  # unregistered input all look like, so the count is printed rather than
  # inferred from finiteness.
  cat(sprintf("\nexactly-zero entries: %d of %d\n",
              sum(g$gradient == 0), length(g$gradient)))
  cat("\ngradient again on the same object (the sweep narrows the system, so a\n")
  cat("second call proves it was widened back)\n")
  g2 <- stand_gradient(scm, traits = smoke_traits)
  cat(sprintf("max |difference| between the two calls: %.3g\n",
              max(abs(g$gradient - g2$gradient))))
}
cat(sprintf("\npeak RSS before the gradient %.1f MB, after %.1f MB\n",
            before_mb, after_mb))

base <- stand_census(scm)
cat("\nre-run central difference\n")
for (tr in smoke_traits) {
  v0 <- TF24_Strategy()$pars[[tr]]
  h <- smoke_fd_eps * abs(v0)
  fd <- (stand_census(smoke_scm(lifetime, tr, v0 + h)) -
         stand_census(smoke_scm(lifetime, tr, v0 - h))) / (2 * h)
  for (m in names(base)) {
    if (inherits(g, "error")) {
      cat(sprintf("%-20s %-6s difference %14.6g\n", m, tr, fd[[m]]))
    } else {
      a <- g$gradient[m, tr]
      cat(sprintf("%-20s %-6s adjoint %14.6g   difference %14.6g   rel %10.3g\n",
                  m, tr, a, fd[[m]], abs(a - fd[[m]]) / max(abs(a), abs(fd[[m]]))))
    }
  }
}
