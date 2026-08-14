## The frontier: consumer error against knot count, on real stands, through
## stand development.
##
## Every earlier comparison fixed a knot count and asked which placement was most
## accurate. That is the wrong question, because the placements do not cost the
## same at the same count -- a canopy-tied grid holds its count fixed while a
## lattice grows one. The question a design has to answer is the frontier: at a
## given cost, which placement is most accurate, at every stage of a run.
##
## The stands are real: states taken from an FF16 run rather than constructed, so
## the height distribution is whatever self-thinning actually produces.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(240)

ff16_scm <- function(lifetime) {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- lifetime
  p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"),
                       hyperpar = FF16_hyperpar, birth_rate = list(20))
  scm <- SCM("FF16", "FF16_Env")(p1, Environment("FF16"), Control())
  scm$collect <- FALSE
  scm$run()
  ## A species does not expose its strategy, so the crown shape comes from the
  ## parameters the run was built from.
  list(scm = scm, eta = p1$strategies[[1]]$pars[["eta"]])
}
heights_of <- function(p) {
  h <- unlist(lapply(p$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))
  sort(h[is.finite(h) & h > 0])
}

## ---- placements, all as functions of the canopy top -------------------------
## A lattice whose spacing grows with height: z_j = z0 r^j above a flat foot.
## Every position is a constant of the run, exactly as a uniform lattice's are,
## so the dropped position channel is zero for the same reason -- but the
## resolution a cohort gets near its own top no longer depends on how tall the
## stand has become.
## z_j = z_c (r^j - 1). This starts at exactly 0, has spacing (z + z_c)(r - 1)
## -- proportional to height well above the foot and constant well below it --
## and inverts in closed form, j = log(1 + z/z_c)/log r, so a query still finds
## its span by arithmetic rather than by a search.
place_graded <- function(hmax, r, zc = 0.3, pad = 2) {
  J <- ceiling(log(1 + hmax / zc) / log(r)) + pad
  zc * (r^(0:J) - 1)
}
place_dyadic <- function(hmax, dmin = 0.0125, target = 200) {
  m <- 0
  while (hmax / (dmin * 2^m) + 2 > target) m <- m + 1
  seq(0, by = dmin * 2^m, length.out = ceiling(hmax / (dmin * 2^m)) + 2)
}

PLACES <- list(
  `develop adaptive`   = function(p, hm, hs) place_adaptive(p, hm),
  `canopy-tied 65`     = function(p, hm, hs) place_uniform_hmax(hm, 65),
  `canopy-tied 129`    = function(p, hm, hs) place_uniform_hmax(hm, 129),
  `canopy-tied 257`    = function(p, hm, hs) place_uniform_hmax(hm, 257),
  `uniform 0.100`      = function(p, hm, hs) place_fixed_abs(hm, 0.100),
  `uniform 0.050`      = function(p, hm, hs) place_fixed_abs(hm, 0.050),
  `uniform 0.025`      = function(p, hm, hs) place_fixed_abs(hm, 0.025),
  `graded r=1.04`      = function(p, hm, hs) place_graded(hm, 1.04),
  `graded r=1.02`      = function(p, hm, hs) place_graded(hm, 1.02),
  `graded r=1.01`      = function(p, hm, hs) place_graded(hm, 1.01),
  `dyadic <=200`       = function(p, hm, hs) place_dyadic(hm)
)

## Relative error of the crown-mean light each cohort actually reads.
consumer_err <- function(patch, kx, hs, eta) {
  ex <- function(q) exact_field(patch, q, "L")$value
  me <- vapply(hs, function(h) crown_mean(ex, h, eta, gl), numeric(1))
  f <- fit_at(patch, kx, "L")
  mf <- vapply(hs, function(h)
    crown_mean(function(q) herm_eval(f, q)$value, h, eta, gl), numeric(1))
  abs(mf - me) / me
}

for (lt in c(2, 4, 6, 10, 20, 40)) {
  run <- ff16_scm(lt)
  p <- run$scm$patch
  hs <- heights_of(p)
  if (!length(hs)) { cat(sprintf("lifetime %g: no cohorts\n", lt)); next }
  hm <- max(hs)
  eta <- run$eta
  stopifnot(length(eta) == 1L, is.finite(eta))
  cat(sprintf("\n=== lifetime %g yr : %d cohorts, canopy %.2f m, median height %.2f m ===\n",
              lt, length(hs), hm, median(hs)))
  cat(sprintf("%-18s %6s %8s %7s | %10s %10s %10s\n", "placement", "knots", "top d/h",
              "min in", "max err", "err tall", "err shortest"))
  for (nm in names(PLACES)) {
    kx <- tryCatch(PLACES[[nm]](p, hm, hs), error = function(e) NULL)
    if (is.null(kx)) { cat(sprintf("%-18s  failed\n", nm)); next }
    e <- consumer_err(p, kx, hs, eta)
    ## spacing of the span holding the canopy top, relative to the canopy
    k <- max(which(kx < hm))
    dtop <- (kx[k + 1] - kx[k]) / hm
    ## the worst-served cohort: how many knots sit in the band where its own
    ## leaf area actually is
    nb <- min(vapply(hs, function(h) sum(kx > 0.7 * h & kx <= h), numeric(1)))
    cat(sprintf("%-18s %6d %8.4f %7d | %10.2e %10.2e %10.2e\n",
                nm, length(kx), dtop, nb, max(e), e[[length(e)]], e[[1]]))
  }
}
