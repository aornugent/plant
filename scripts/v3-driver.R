# V3: one step's lambda_y against a central difference of the same step, one
# seeded output component at a time, over the six strategy-rate columns of every
# node.
#
#   Rscript scripts/v3-driver.R <plant worktree> [sweep]
#
# Pass a second argument to sweep the difference step instead of printing the
# per-row table.
args <- commandArgs(TRUE)
tree <- args[1]
if (is.na(tree)) stop("give the plant worktree as the first argument")
do_sweep <- length(args) > 1

suppressMessages({ library(odelia); pkgload::load_all(tree, quiet = TRUE, export_all = TRUE) })
Sys.setenv(TESTTHAT_PARALLEL = "false")
Sys.setenv(PKG_CPPFLAGS = paste0("-I", tree, "/inst/include -DNDEBUG"),
           PKG_LIBS = file.path(tree, "src/plant.so"))
Rcpp::sourceCpp(file.path(tree, "scratch/wire_gates.cpp"))

# Two control values differ from the defaults, and both belong to the reference
# rather than to the adjoint. Read at the defaults this comparison reports the
# adjoint as four orders too small on the log_density rows, and that reading is
# an artefact of the difference, not a defect of the transpose:
#
#   GSS_tol_abs. The collar-potential bracket search returns a point that moves
#   in steps as the bracket lands differently, so at the default 1e-1 one step's
#   y_end is not Lipschitz at the difference step. Measured at 1e-5
#   displacements of one cohort's storage: y_end log_density moves by 1.14e-03
#   where the true derivative accounts for 9.2e-13, and y_end height moves by
#   4.29e-06. At 1e-6 those become 1.59e-10 and 7.55e-12.
#
#   node_gradient_eps. The transport term is a difference quotient divided by
#   this, so it multiplies whatever irreproducibility the growth rate carries by
#   its reciprocal. At the default 1e-6 that is a factor of a million and it is
#   what puts the log_density rows four orders out.
#
# Both are forward settings, so the adjoint and the difference are taken at the
# same operator and the comparison stays honest; what changes is only whether a
# difference of the forward model can resolve it.
v3_gss_tol_abs <- 1e-6
v3_node_gradient_eps <- 1e-3

# The step size is inside the sizes this run accepts, so it is a step the model
# takes. frac 0.5 leaves the youngest cohort a growth interval: a cohort
# introduced at the end of the run is still at height_0 with a zero-width
# interval below it and its rows dominate any readout over the whole state.
v3_step_size <- 0.02
v3_fd_eps <- 1e-7

v3_patch <- function(lifetime = 20, n_intro = 8, frac = 0.5) {
  ctrl <- Control()
  ctrl$GSS_tol_abs <- v3_gss_tol_abs
  ctrl$node_gradient_eps <- v3_node_gradient_eps
  s <- TF24_Strategy()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  p <- scm_base_parameters("TF24")
  p$strategies <- list(s)
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list(seq(0, lifetime * frac, length.out = n_intro))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), ctrl)
  scm$run()
  scm$patch
}

# Patch::ode_state writes the species first and the environment last, so a node
# occupies eight consecutive slots -- six strategy states, offspring, log
# density -- and the environment slots trail.
node_stride <- 8L
n_strategy_states <- 6L

# The columns the comparison is read over. The difference returns a column for
# every state including offspring, log density and the environment, and those
# are left out here for the reason V1 leaves them out.
strategy_columns <- function(n, nodes) {
  i0 <- seq_len(n) - 1L
  which(i0 < node_stride * nodes & (i0 %% node_stride) < n_strategy_states)
}

# Both figures are reported for every row because they diagnose differently.
# Pointwise is |a - b| / max(|a|, |b|) entry by entry, and an entry whose true
# value sits below the reference's own resolution reads 1 whatever the adjoint
# does. Normwise divides the largest entry disagreement by the largest entry of
# the reference, so it asks whether the row is right at the row's own scale.
pointwise <- function(a, b) {
  s <- pmax(abs(a), abs(b))
  max(ifelse(s == 0, 0, abs(a - b) / s))
}
normwise <- function(a, b) max(abs(a - b)) / max(abs(b))

fd_jacobian <- function(patch, y, t0, eps) {
  n <- length(y)
  sapply(seq_len(n), function(j) {
    up <- y; up[j] <- up[j] + eps
    down <- y; down[j] <- down[j] - eps
    (step_forward(patch, t0, v3_step_size, up) -
     step_forward(patch, t0, v3_step_size, down)) / (2 * eps)
  })
}

adjoint_rows <- function(patch, y, t0, rows) {
  n <- length(y)
  t(sapply(rows, function(i) {
    seed <- numeric(n); seed[i] <- 1
    step_adjoint_gate(patch, t0, v3_step_size, y, seed)$lambda_in
  }))
}

patch <- v3_patch()
n <- patch$ode_size
nodes <- length(patch$species[[1]]$nodes)
heights <- sapply(patch$species[[1]]$nodes,
                  function(z) z$individual$state("height"))
y <- patch$ode_state
t0 <- patch$time
cols <- strategy_columns(n, nodes)
rows <- seq_len(node_stride * nodes)
adj <- adjoint_rows(patch, y, t0, rows)

cat(sprintf("patch: time %.4g  ode_size %d  nodes %d  step %.4g\n",
            t0, n, nodes, v3_step_size))
cat("heights:", sprintf("%.4f", heights), "\n")

if (do_sweep) {
  transport <- (seq_len(nodes) - 1L) * node_stride + node_stride
  cat(sprintf("%10s %14s %14s %14s %14s\n", "fd_eps", "normwise", "  of which",
              "  the rest", "pointwise"))
  cat(sprintf("%10s %14s %14s %14s %14s\n", "", "max", "transport", "", "max"))
  for (eps in 10^seq(-8, -5, by = 0.5)) {
    fd <- fd_jacobian(patch, y, t0, eps)
    nw <- sapply(rows, function(i) normwise(adj[i, cols], fd[i, cols]))
    pw <- sapply(rows, function(i) pointwise(adj[i, cols], fd[i, cols]))
    cat(sprintf("%10.2e %14.3g %14.3g %14.3g %14.3g\n", eps, max(nw),
                max(nw[transport]), max(nw[-transport]), max(pw)))
  }
} else {
  fd <- fd_jacobian(patch, y, t0, v3_fd_eps)
  cat(sprintf("difference step %.4g\n\n", v3_fd_eps))
  cat(sprintf("%-14s %12s %12s %12s %12s\n", "seeded row", "normwise",
              "pointwise", "|adjoint|", "|difference|"))
  for (k in seq_len(nodes)) {
    for (sl in seq_len(node_stride)) {
      i <- (k - 1L) * node_stride + sl
      cat(sprintf("node %d slot %d  %12.3g %12.3g %12.4g %12.4g\n", k, sl,
                  normwise(adj[i, cols], fd[i, cols]),
                  pointwise(adj[i, cols], fd[i, cols]),
                  sqrt(sum(adj[i, cols]^2)), sqrt(sum(fd[i, cols]^2))))
    }
  }
  nw <- sapply(rows, function(i) normwise(adj[i, cols], fd[i, cols]))
  cat(sprintf("\nnormwise over all %d rows: max %.3g at row %d\n",
              length(rows), max(nw), which.max(nw)))
}

# Measured on plant 4fff1e22 with the odelia holding Step::step_adjoint over
# AdjointRates. At the difference step above, every one of the 64 rows agrees
# normwise, worst 1.36e-02 at node 2 slot 8; the plateau runs from 1e-08 to
# 3.16e-07 and the reading degrades to 0.104 at 1e-06 and 1.3 at 1e-05.
# Pointwise the same table reads up to 1.87, on entries the difference cannot
# resolve: its floor over the rows that carry no transport term is 1.26e-03 of
# the row maximum, and those entries are orders below it.
