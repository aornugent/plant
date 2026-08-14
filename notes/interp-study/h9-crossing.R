## Two things coincide when the canopy top sits on a knot: the knot COUNT
## changes (ceil straddles), and the tallest cohort's curvature BREAK crosses a
## knot. They are separable by holding the count fixed.
##
## The quantity swept is the crown-mean light the tallest cohort reads, as a
## function of its own height, which is what its physiology consumes.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

gl <- gl_nodes(220)
D  <- 0.25

stand_at <- function(h_top) {
  hs <- c(h_top, seq(15.5, 0.6, length.out = 7))
  ld <- -2.2 + seq(0, 0.7, length.out = 8)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

## grid built by the ceil rule: the count moves as the canopy grows
grid_ceil  <- function(h) seq(0, by = D, length.out = ceiling(h / D) + 3)
## grid of a fixed count, generous enough to cover the whole sweep: the count
## cannot move, so only the break crossing is left
grid_fixed <- function(h) seq(0, by = D, length.out = 100)

sweep <- function(grid_fn, label) {
  ## walk h_top across the knot at 18.25
  hs <- 18.25 + seq(-0.02, 0.02, length.out = 21)
  m <- vapply(hs, function(h) {
    p <- stand_at(h)
    eta <- ladder_strategy_parameter(p, 1, "eta")
    f <- fit_at(p, grid_fn(h), "L")
    ff <- function(q) herm_eval(f, q)$value
    ex <- function(q) exact_field(p, q, "L")$value
    c(fit = crown_mean(ff, h, eta, gl), exact = crown_mean(ex, h, eta, gl),
      nk = length(grid_fn(h)))
  }, numeric(3))
  d_fit <- diff(m["fit", ]) / diff(hs)
  d_ex  <- diff(m["exact", ]) / diff(hs)
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("knot count over the sweep : %d .. %d\n",
              min(m["nk", ]), max(m["nk", ])))
  cat(sprintf("crown mean, max |fit-exact|: %.3e\n",
              max(abs(m["fit", ] - m["exact", ]))))
  ## a kink shows as a jump in the slope of the error, not of the value
  d_err <- d_fit - d_ex
  cat(sprintf("d(error)/dh  range          : [%.3e, %.3e]  jump %.3e\n",
              min(d_err), max(d_err), max(abs(diff(d_err)))))
  cat(sprintf("value continuous across knot: max |second difference| %.3e\n",
              max(abs(diff(diff(m["fit", ]))))))
  invisible(NULL)
}

sweep(grid_ceil,  "ceil rule: count changes at the crossing")
sweep(grid_fixed, "fixed count: only the break crosses")

## And for contrast, the canopy-following grid, where the top knot IS h_top
## always, so the break never crosses -- but the positions all move.
sweep(function(h) seq(0, h, length.out = 65),
      "canopy-following: break pinned to the top knot, positions move")
