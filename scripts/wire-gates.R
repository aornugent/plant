# Build and run the gate harness in scratch/wire_gates.cpp, at states pinned
# here. Run it with the plant source tree already built:
#
#   R_LIBS_USER=<a library holding odelia> Rscript scripts/wire-gates.R
#
# A gate whose state lives in a transcript is not re-runnable. Everything the
# readouts below depend on -- the compile flags, the patch, the seed -- is in
# this file.

plant_dir <- Sys.getenv("PLANT_DIR", unset = normalizePath("."))

suppressMessages(pkgload::load_all(plant_dir, quiet = TRUE, export_all = TRUE))
suppressMessages(library(odelia))
Sys.setenv(TESTTHAT_PARALLEL = "false")

# 1. The compile.
#
# Two things beyond the usual AD-test recipe. The snippet calls into plant's
# own templates, which are not all header-inline, so PKG_LIBS needs plant.so as
# well as odelia.so -- without it the load fails on
# plant::quadrature::QK::integrate_vector. And the driver reaches Parameters,
# Patch and SCM, which load_all only exports under export_all = TRUE.
wire_gates_compile <- function(dir = plant_dir) {
  plant_so <- file.path(dir, "src/plant.so")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  stopifnot(file.exists(plant_so), file.exists(odelia_so))
  Sys.setenv(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(file.path(dir, "inst/include"))),
                         paste0("-I", shQuote(system.file("include", package = "odelia")))),
    PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                     shQuote(normalizePath(odelia_so))))
  Rcpp::sourceCpp(file.path(dir, "scratch/wire_gates.cpp"),
                  env = globalenv(), verbose = FALSE)
}

# 2. The state.
#
# scm_base_parameters("TF24") returns strategies = list(), so a strategy has to
# be assigned before the Patch constructor sees it; taking the pattern from
# test-patch.R instead gives "subscript out of bounds". The run is short and its
# schedule explicit, so the same call gives the same state every time.
#
# frac places the last introduction. A cohort introduced at the end of the run
# is still at height_0 with the interval below it of zero width, which is a
# degenerate configuration rather than a neutral one: its rows dominate any
# readout taken over the whole state. frac = 0.5 leaves the youngest cohort a
# growth interval.
harness_scm <- function(lifetime = 20, n_intro = 8, frac = 0.5) {
  s <- TF24_Strategy()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- scm_base_parameters("TF24")
  p$strategies <- list(s)
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list(seq(0, lifetime * frac, length.out = n_intro))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  scm$run()
  scm
}

# 3. The state layout, and the seed.
#
# Species first, environment last: Patch::ode_state walks the nodes and then the
# environment, so with a zero-based index
#
#   node        = index / node_stride
#   slot        = index % node_stride
#   environment iff index >= node_stride * node_count
#
# and there is no leading offset. A node's eight slots are its six strategy
# states, then offspring, then log density.
#
# The readout seeds only the six strategy-rate slots, and compares only the six
# strategy columns. Both halves of that are needed. Seeding a rate the
# decomposition does not carry -- log_density carries the transport term,
# offspring the fecundity weighting -- puts a channel in the reference that is
# not in the accumulator. And reading the residual over all columns pins it at
# exactly 1 whatever the state: seeding a subset of rows does not remove
# columns, so the whole-patch recording returns nonzero adjoints in the
# log_density and environment columns where the decomposition returns zero by
# construction. Restricted to the columns both models contain, the readout
# discriminates.
node_stride <- 8L
n_strategy_states <- 6L

node_slot <- function(k, s) (k - 1L) * node_stride + s
transport_slot <- function(k) node_slot(k, node_stride)

rate_only_seed <- function(patch, cohorts = NULL) {
  n_node <- length(patch$species[[1]]$nodes)
  lam <- numeric(patch$ode_size)
  if (is.null(cohorts)) cohorts <- seq_len(n_node)
  for (k in cohorts) {
    lam[node_slot(k, seq_len(n_strategy_states))] <-
      round(seq(0.3, 1.7, length.out = n_strategy_states) * c(1, -1), 6)
  }
  lam
}

# The columns both the decomposition and the whole-patch recording express.
strategy_columns <- function(patch) {
  n_node <- length(patch$species[[1]]$nodes)
  as.vector(sapply(seq_len(n_node), function(k) node_slot(k, seq_len(n_strategy_states))))
}

rel <- function(a, b) sqrt(sum((a - b)^2)) / sqrt(sum(b^2))

# 4. The V3 state and seed: one step from the pinned patch, one seeded component
# at a time. The step size is well inside the accepted sizes of the run above,
# so the step is one the model actually takes.
v3_step_size <- 0.02
v3_fd_eps <- 1e-6

v1_readout <- function(patch, cohorts = NULL) {
  lam <- rate_only_seed(patch, cohorts)
  cols <- strategy_columns(patch)
  ref <- patch_recording_vjp(patch, lam)[cols]
  sapply(1:5, function(u)
    rel(patch_adjoint_partial(patch, lam, u)$lambda_y[cols], ref))
}

# The transport row of one cohort's block, against a central difference of the
# same block. The reference must be a difference of the quotient the model
# solves, never an analytic dg/dh: the sub-grid probe and the derivative it
# approximates are different operators here. The agreement is reported beside
# the absolute scale, because the quotient differences two nearly-equal
# parameter derivatives and divides by node_gradient_eps, which leaves about
# 1e-10 of absolute error before any question of smoothness.
transport_block_gate <- function(patch, species_index = 0, node_index = 3,
                                 inputs = c(1, 2, 5, 7, 20, 90)) {
  b <- block_vjp(patch, species_index, node_index, rep(0, 12))
  x0 <- b$inputs
  row <- n_strategy_states + 1L          # the transport output
  seed <- rep(0, b$n_out); seed[row] <- 1
  adj <- block_vjp(patch, species_index, node_index, seed)$adjoint
  inputs <- inputs[inputs <= length(x0)]
  t(sapply(inputs, function(j) {
    h <- max(1e-7, abs(x0[j]) * 1e-7)
    up <- x0; up[j] <- up[j] + h
    down <- x0; down[j] <- down[j] - h
    d <- (block_value(patch, species_index, node_index, up)[row] -
          block_value(patch, species_index, node_index, down)[row]) / (2 * h)
    c(input = j, adjoint = adj[j], difference = d, abs_err = abs(adj[j] - d))
  }))
}

if (sys.nframe() == 0L) {
  wire_gates_compile()
  patch <- harness_scm()$patch
  n_node <- length(patch$species[[1]]$nodes)
  cols <- strategy_columns(patch)

  cat("state: time", patch$time, " ode_size", patch$ode_size,
      " cohorts", n_node, "\n")
  cat("patch reaches the sweep through its own rate transpose:",
      patch_carries_rates_adjoint(), "\n\n")

  cat("V1-shaped readout on the strategy columns, one contribution at a time\n")
  cat("  (a)soil (a)+offspring (b)+blocks (c)+knots (d)+allometry\n   ",
      sprintf("%.4g", v1_readout(patch)), "\n\n")

  cat("transport row of the block, against a difference of the same quotient\n")
  print(signif(transport_block_gate(patch), 6))

  lam <- rate_only_seed(patch)
  lam_t <- lam
  for (k in seq_len(n_node)) lam_t[transport_slot(k)] <- 0.5
  with_t <- patch_adjoint(patch, lam_t)$lambda_y
  # Severed: the transport rate adjoint dropped rather than seeded, which is
  # what the patch did before the block carried the term.
  without_t <- patch_adjoint(patch, lam)$lambda_y
  cat("\ntransport contribution to lambda_y: max|.| =",
      sprintf("%.6g", max(abs(with_t - without_t))),
      " severed:", sprintf("%.6g", max(abs(without_t - without_t))), "\n")

  cat("\nV3: one step's lambda_y on the strategy columns, against a finite\n")
  cat("difference of the same step, one seeded rate component at a time\n")
  y <- patch$ode_state
  t0 <- patch$time
  fd <- sapply(seq_along(y), function(j) {
    up <- y; up[j] <- up[j] + v3_fd_eps
    down <- y; down[j] <- down[j] - v3_fd_eps
    (step_forward(patch, t0, v3_step_size, up) -
     step_forward(patch, t0, v3_step_size, down)) / (2 * v3_fd_eps)
  })
  for (k in seq_len(n_node)) {
    for (sl in seq_len(node_stride)) {
      i <- node_slot(k, sl)
      seed <- numeric(length(y)); seed[i] <- 1
      got <- step_adjoint_gate(patch, t0, v3_step_size, y, seed)$lambda_in
      cat(sprintf("  node %d slot %d  rel %.4g\n", k, sl,
                  rel(got[cols], fd[i, cols])))
    }
  }
}
