## Interpolate the optical depth, or the light?
##
## The shipped field applies Beer's law BEFORE interpolation: the knots hold
## L = exp(-A) and dL/dz = -A' exp(-A). The alternative holds A and A' and
## exponentiates at the query.
##
## This is not a placement question and it is not cosmetic. What a cohort
## consumes is a RELATIVE error in the light it reads, and relative error in
## L = exp(-A) is absolute error in A. In deep shade L spans decades while A
## stays O(1), so a knot budget that holds L to a fixed absolute accuracy holds
## the suppressed cohorts -- the ones self-thinning kills -- to a far worse
## relative one. A is also the quantity the reduction actually forms, so holding
## it removes an exponential from the build and moves it to the query.
##
## Both halves matter, so both are measured: the spacing each cohort requires,
## and the count a placement needs to meet it.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(240)
TARGET <- 1e-6

ff16_scm <- function(lifetime) {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- lifetime
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                       hyperpar = FF16_hyperpar, birth_rate = list(20))
  scm <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
  scm$collect <- FALSE
  scm$run()
  list(patch = scm$patch, eta = p1$strategies[[1]]$pars[["eta"]])
}

## The consumer is the same in both spaces: the crown-mean of the LIGHT. Only
## what the knots hold differs, so the comparison is like for like.
light_of <- function(f, space) {
  if (space == "L") function(q) herm_eval(f, q)$value
  else function(q) exp(-herm_eval(f, q)$value)
}

err_at <- function(patch, h, eta, kx, space) {
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- crown_mean(ex, h, eta, gl)
  f <- fit_at(patch, kx, space)
  mf <- crown_mean(light_of(f, space), h, eta, gl)
  abs(mf - me) / me
}

required <- function(patch, h, eta, hmax, space) {
  grid <- function(d) seq(0, by = d, length.out = ceiling(hmax / d) + 3)
  lo <- 1e-3; hi <- 0.5
  if (err_at(patch, h, eta, grid(lo), space) > TARGET) return(NA_real_)
  if (err_at(patch, h, eta, grid(hi), space) < TARGET) return(hi)
  for (i in 1:18) {
    mid <- sqrt(lo * hi)
    if (err_at(patch, h, eta, grid(mid), space) < TARGET) lo <- mid else hi <- mid
  }
  lo
}

for (lt in c(4, 10, 40)) {
  run <- ff16_scm(lt)
  p <- run$patch; eta <- run$eta
  hs <- sort(unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1)))))
  hs <- hs[is.finite(hs) & hs > 0]
  hmax <- max(hs)
  Lg <- exact_field(p, c(0, hmax / 2), "L")$value
  cat(sprintf("\n=== lifetime %g yr, canopy %.2f m, L(0)=%.3g ===\n", lt, hmax, Lg[[1]]))

  cat("\n  spacing each cohort requires, by what the knots hold:\n")
  cat(sprintf("  %8s %12s %12s %8s\n", "height", "holding L", "holding A", "gain"))
  idx <- unique(round(seq(1, length(hs), length.out = 8)))
  for (i in idx) {
    h <- hs[[i]]
    dL <- required(p, h, eta, hmax, "L")
    dA <- required(p, h, eta, hmax, "A")
    cat(sprintf("  %8.3f %12.4f %12.4f %7.2fx\n", h, dL, dA, dA / dL))
  }

  cat("\n  consumer error at a fixed uniform spacing, both spaces:\n")
  cat(sprintf("  %8s %7s | %11s %11s | %11s %11s\n", "spacing", "knots",
              "L max", "L shortest", "A max", "A shortest"))
  for (d in c(0.1, 0.05, 0.025)) {
    kx <- seq(0, by = d, length.out = ceiling(hmax / d) + 3)
    eL <- vapply(hs, function(h) err_at(p, h, eta, kx, "L"), numeric(1))
    eA <- vapply(hs, function(h) err_at(p, h, eta, kx, "A"), numeric(1))
    cat(sprintf("  %8.3f %7d | %11.2e %11.2e | %11.2e %11.2e\n",
                d, length(kx), max(eL), eL[[1]], max(eA), eA[[1]]))
  }
}
