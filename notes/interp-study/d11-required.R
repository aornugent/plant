## What spacing does each cohort actually need?
##
## The frontier run says a uniform lattice starves the short cohorts and a
## geometric one starves the tall, so neither shape is right and picking between
## them is picking which end to fail at. Rather than argue a third shape, measure
## the thing a shape has to match: for a cohort of height h, the knot spacing at
## which the crown-mean light it reads is accurate to a stated target.
##
## Delta_required(h) IS the design curve. A fixed grid is good exactly insofar as
## its spacing sits under that curve everywhere, and its cost is the integral of
## 1/Delta over the column. Nothing here is a placement -- it is the requirement
## every placement is then judged against.
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

## One cohort's crown-mean error under a uniform lattice of spacing d. Uniform,
## because the question is what spacing this cohort needs LOCALLY -- a global
## shape is what gets designed afterwards from the answer.
err_at <- function(patch, h, eta, d, hmax) {
  kx <- seq(0, by = d, length.out = ceiling(hmax / d) + 3)
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- crown_mean(ex, h, eta, gl)
  f <- fit_at(patch, kx, "L")
  mf <- crown_mean(function(q) herm_eval(f, q)$value, h, eta, gl)
  abs(mf - me) / me
}

## Bisect on log spacing for the coarsest d meeting the target.
required <- function(patch, h, eta, hmax) {
  lo <- 1e-3; hi <- 0.5
  if (err_at(patch, h, eta, lo, hmax) > TARGET) return(NA_real_)
  if (err_at(patch, h, eta, hi, hmax) < TARGET) return(hi)
  for (i in 1:18) {
    mid <- sqrt(lo * hi)
    if (err_at(patch, h, eta, mid, hmax) < TARGET) lo <- mid else hi <- mid
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
  ## a spread of cohorts, not all 100 of them
  idx <- unique(round(seq(1, length(hs), length.out = 14)))
  cat(sprintf("\n=== lifetime %g yr, canopy %.2f m, target %.0e ===\n",
              lt, hmax, TARGET))
  cat(sprintf("%8s %12s %12s %10s\n", "height", "d required", "d/h", "knots if uniform"))
  for (i in idx) {
    h <- hs[[i]]
    d <- required(p, h, eta, hmax)
    cat(sprintf("%8.3f %12.4f %12.5f %10s\n", h, d, d / h,
                if (is.na(d)) "-" else sprintf("%.0f", hmax / d)))
  }
}
