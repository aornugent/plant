## Does the reduction really have the shape the sweep needs?
##
## The sweep rests on the claim that the field is a weighted sum over cohorts,
##
##   A(z) = sum_k c_k Q(z / h_k),   c_k = w_k n_k k_I a_k
##
## with w_k the trapezium weights on the birth-date abscissa and c_k independent
## of z. If that is wrong the whole design is wrong, and it is cheaper to find
## out here than in C++.
##
## Two checks, in order. First rebuild the model's own reduction from per-node
## quantities and compare -- that tests the claim about the weights. Then run
## the sweep against the rebuild -- that tests the expansion.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
suppressMessages({
  library(odelia)
  pkgload::load_all("/home/a/dev/plant-dev/plant/.claude/worktrees/fixed-grid", quiet = TRUE)
})

p0 <- scm_base_parameters("FF16"); p0$max_patch_lifetime <- 10
p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                     hyperpar = FF16_hyperpar, birth_rate = list(20))
s <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
s$collect <- FALSE; s$run()
p <- s$patch; sp <- p$species[[1]]

eta <- p1$strategies[[1]]$pars[["eta"]]
k_I <- p1$strategies[[1]]$pars[["k_I"]]
x <- sp$node_times                       # the birth-date abscissa
h <- sp$heights
a <- sp$compute_competition_effect_by_nodes
n <- exp(sp$log_densities)
cat(sprintf("species: %d nodes, abscissa %.4f..%.4f, heights %.3f..%.3f\n",
            length(x), min(x), max(x), min(h), max(h)))

## Trapezium weights on the abscissa, in the order the nodes are stored. The
## model accumulates (x0 - x1)(f1 + f0) / 2 over consecutive nodes, which is
## this sum re-associated.
trap_weights <- function(x) {
  m <- length(x)
  if (m == 1) return(0)
  w <- numeric(m)
  d <- diff(x)
  w[1] <- d[1] / 2
  w[m] <- d[m - 1] / 2
  if (m > 2) w[2:(m - 1)] <- (x[3:m] - x[1:(m - 2)]) / 2
  w
}
w <- trap_weights(x)
cc <- w * n * k_I * a                     # the per-cohort amplitude

Q <- function(u) ifelse(u > 1, 0, (1 - u^eta)^2)
dQdz <- function(z, hh) ifelse(z > hh, 0, (-2 * eta * (z / hh)^eta +
                                            2 * eta * (z / hh)^(2 * eta)) / z)

rebuild <- function(z) {
  v <- vapply(z, function(zz) sum(cc * Q(zz / h)), numeric(1))
  sl <- vapply(z, function(zz) if (zz <= 0) 0 else sum(cc * dQdz(zz, h)), numeric(1))
  list(value = v, slope = sl)
}

## The model's own reduction, boundary interval left off so it is the sum over
## the stored nodes alone -- which is what the rebuild covers.
zz <- seq(0.05, max(h) * 0.999, length.out = 60)
model <- vapply(zz, function(q) p$compute_competition_and_slope(q), numeric(2))
rb <- rebuild(zz)
cat(sprintf("\nrebuild against the model's reduction (boundary interval included in the model):\n"))
cat(sprintf("  value: max abs diff %.3e over a range of %.3f\n",
            max(abs(rb$value - model[1, ])), diff(range(model[1, ]))))
cat(sprintf("  slope: max abs diff %.3e over a range of %.3f\n",
            max(abs(rb$slope - model[2, ])), diff(range(model[2, ]))))
cat("  (a residual of the boundary node's own trapezium is expected here)\n")

## The sweep, against the rebuild. Same weights, same amplitudes, so any
## difference is the expansion and its conditioning alone.
sweep_field <- function(z) {
  ord <- order(z, decreasing = TRUE)
  zs <- z[ord]
  oh <- order(h, decreasing = TRUE)
  hs <- h[oh]; cs <- cc[oh]
  v <- numeric(length(zs)); sl <- numeric(length(zs))
  S0 <- 0; T1 <- 0; T2 <- 0; prev <- NA_real_; j <- 1L
  for (i in seq_along(zs)) {
    zi <- zs[i]
    if (!is.na(prev)) { r <- (zi / prev)^eta; T1 <- T1 * r; T2 <- T2 * r * r }
    while (j <= length(hs) && hs[j] >= zi) {
      u <- (zi / hs[j])^eta
      S0 <- S0 + cs[j]; T1 <- T1 + cs[j] * u; T2 <- T2 + cs[j] * u * u
      j <- j + 1L
    }
    v[i] <- S0 - 2 * T1 + T2
    sl[i] <- if (zi > 0) (-2 * eta * T1 + 2 * eta * T2) / zi else 0
    prev <- zi
  }
  o <- integer(length(zs)); o[ord] <- seq_along(zs)
  list(value = v[o], slope = sl[o])
}

zk <- seq(0, max(h) * 1.02, by = 0.05)    # the lattice this branch builds on
rb2 <- rebuild(zk); sw <- sweep_field(zk)
cat(sprintf("\nsweep against the rebuild, on the shipped lattice (%d knots):\n", length(zk)))
cat(sprintf("  A(0) = %.4f, so the direct sum's own rounding is about %.2e\n",
            rb2$value[[1]], .Machine$double.eps * rb2$value[[1]]))
cat(sprintf("  value: max abs diff %.3e\n", max(abs(rb2$value - sw$value))))
cat(sprintf("  slope: max abs diff %.3e   (slope range %.3f)\n",
            max(abs(rb2$slope - sw$slope)), diff(range(rb2$slope))))
cat(sprintf("\n  kernel evaluations: per-height %d, swept %d, ratio %.0fx\n",
            sum(vapply(zk, function(q) sum(h >= q), numeric(1))),
            length(zk) + length(h),
            sum(vapply(zk, function(q) sum(h >= q), numeric(1))) /
              (length(zk) + length(h))))
