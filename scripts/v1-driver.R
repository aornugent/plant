# V1: the hand-written reverse-pass decomposition against one whole-Patch
# recording, added one contribution at a time. Builds the patch at an explicit
# state, seeds the adjoint, and prints the five-row incremental readout.
#
# Rscript scripts/v1-driver.R [max_patch_lifetime]
library(odelia)
pkgload::load_all("/home/user/wt-gt", quiet = TRUE, export_all = TRUE)

Sys.setenv(PKG_CPPFLAGS = "-I/home/user/wt-gt/inst/include -DNDEBUG",
           PKG_LIBS = "/home/user/wt-gt/src/plant.so")
Rcpp::sourceCpp("/home/user/wt-gt/scratch/wire_gates.cpp")

# State layout of the patch ode vector: node_count blocks of 8, then 9 trailing
# environment slots. Within a node: slots 0-5 the strategy states, 6 offspring,
# 7 log_density. (Patch::ode_state writes the species first, environment last.)
v1_patch <- function(lifetime = 2, lma = 0.1978791) {
  p <- scm_base_parameters("TF24", "TF24_Env")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(lma, "lma"))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  scm$run()
  scm$patch
}

# Only the six strategy-rate slots of each node are seeded. Two channels the
# recording carries and the decomposition does not are excluded here rather
# than toleranced later: the soil (the soil store is a declared passive
# boundary, so the decomposition has no soil channel at all) and the transport
# stencil (reached through log_density, slot 7, omitted until P3.5).
v1_seed <- function(patch) {
  n <- patch$ode_size
  nodes <- patch$species[[1]]$size
  full <- as.numeric(seq_len(n)) / n
  lam <- rep(0, n)
  for (k in seq_len(nodes)) for (s in 1:6) {
    i <- (k - 1) * 8 + s
    lam[i] <- full[i]
  }
  lam
}

# patch_adjoint_partial(patch, lambda, upto) runs the decomposition with the
# first `upto` of the five contributions applied, in the order
# ode_rates_adjoint applies them: 1 soil, 2 offspring, 3 cohort blocks,
# 4 light-knot pullback, 5 allometry. It builds its own block seeds from
# lambda and returns lambda_y over the whole ode vector.
v1_readout <- function(patch) {
  n <- patch$ode_size
  nodes <- patch$species[[1]]$size
  lam <- v1_seed(patch)
  ref <- patch_recording_vjp(patch, lam)

  # Compare on the columns both models contain: the six strategy-rate slots.
  # The log_density and environment columns are the two excluded channels and
  # read rel 1 by construction.
  i0 <- seq_len(n) - 1L
  keep <- i0 < 8 * nodes & (i0 %% 8L) <= 5L

  labels <- c("(a) soil", "(a) + offspring", "(b) + cohort blocks",
              "(c) + knot pullback", "(d) + allometry")
  cat(sprintf("patch: time %.4g  ode_size %d  nodes %d\n", patch$time, n, nodes))
  cat(sprintf("%-22s %12s %12s\n", "contributions", "pointwise", "normwise"))
  for (k in 1:5) {
    d <- patch_adjoint_partial(patch, lam, k)$lambda_y
    s <- pmax(abs(d), abs(ref))
    pt <- max(ifelse(s == 0, 0, abs(d - ref) / s)[keep])
    nw <- max(abs((d - ref)[keep])) / max(abs(ref[keep]))
    cat(sprintf("%-22s %12.3g %12.3g\n", labels[k], pt, nw))
  }
}

args <- commandArgs(trailingOnly = TRUE)
lifetime <- if (length(args)) as.numeric(args[[1]]) else 2
v1_readout(v1_patch(lifetime))

# Measured at lifetime 2, lma 0.1978791, seed seq_len(n)/n:
#   (a) soil              1          1
#   (a) + offspring       1          1
#   (b) + cohort blocks   1          1.82e-04
#   (c) + knot pullback   1          1.82e-04
#   (d) + allometry       2.05e-11   3.33e-15
# Pointwise, the readout is not known to close at 1e-11 or better at any state
# tried; 2.05e-11 at lifetime 2 is the best found. Normwise it closes at
# 3.33e-15 here and stays below 1e-11 out to lifetime 3.
