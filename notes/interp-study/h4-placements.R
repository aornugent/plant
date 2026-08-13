## The placement study: accuracy, the dropped position channel, and the
## consumer's own error, for every candidate at a reported knot count.
source("/home/a/.claude/jobs/e02c60e6/tmp/lib-field.R")

shaded_patch <- function(n = 8, top = 18, floor_h = 0.6, ldens = -1.6) {
  hs <- seq(top, floor_h, length.out = n)
  ld <- ldens + seq(0, 0.7, length.out = n)
  ladder_patch(species = "fast", heights = list(hs), log_densities = list(ld))
}

heights_of <- function(patch)
  unlist(lapply(patch$species, function(s)
    vapply(s$nodes, function(n) n$height, numeric(1))))

gl <- gl_nodes(220)

study <- function(patch, label, space = "L") {
  hs   <- heights_of(patch)
  hmax <- max(hs)
  eta  <- ladder_strategy_parameter(patch, 1, "eta")
  z    <- test_grid(hmax, hs)
  ex   <- exact_field(patch, z, space)
  rng  <- diff(range(ex$value))

  ## the consumer: crown-mean light, refereed against the exact field
  exact_fn <- function(zz) exact_field(patch, zz, "L")$value
  mean_exact <- vapply(hs, function(h) crown_mean(exact_fn, h, eta, gl), numeric(1))

  cands <- list(
    `uniform x h_max, 65`   = function(hm) place_uniform_hmax(hm, 65),
    `uniform x h_max, 129`  = function(hm) place_uniform_hmax(hm, 129),
    `uniform x h_max, 257`  = function(hm) place_uniform_hmax(hm, 257),
    `breaks + fill, 65`     = function(hm) place_breaks_fill(hm, hs, 65),
    `fixed abs, d=0.25`     = function(hm) place_fixed_abs(hm, 0.25),
    `fixed abs, d=0.10`     = function(hm) place_fixed_abs(hm, 0.10),
    `fixed geom, r=1.20`    = function(hm) place_fixed_geom(hm, 1.20, 0.02),
    `fixed geom, r=1.10`    = function(hm) place_fixed_geom(hm, 1.10, 0.02),
    `fixed geom, r=1.05`    = function(hm) place_fixed_geom(hm, 1.05, 0.02)
  )
  ad <- place_adaptive(patch, hmax, atol = 1e-6, rtol = 1e-4, space = space)
  cands[[sprintf("adaptive on value, %d", length(ad))]] <- function(hm) {
    place_adaptive(patch, hm, atol = 1e-6, rtol = 1e-4, space = space)
  }

  cat(sprintf("\n=== %s  [space = %s] ===\n", label, space))
  cat(sprintf("h_max %.3f   cohorts %d   L(0) %.4g   field range %.4g\n",
              hmax, length(hs), exact_field(patch, 0, "L")$value, rng))
  cat(sprintf("%-24s %5s %10s %10s %10s %10s %10s %11s\n",
              "placement", "knots", "val max", "val rms", "slope max",
              "slope rms", "min fit", "crown max"))
  for (nm in names(cands)) {
    place <- cands[[nm]]
    kx <- place(hmax)
    f  <- fit_at(patch, kx, space)
    fv <- herm_eval(f, z)
    ev <- err_stats(fv$value, ex$value)
    es <- err_stats(fv$slope, ex$slope)
    ## consumer error: the crown mean each cohort would actually read
    fit_fn <- if (space == "A") function(zz) exp(-herm_eval(f, zz)$value)
              else              function(zz) herm_eval(f, zz)$value
    mean_fit <- vapply(hs, function(h) crown_mean(fit_fn, h, eta, gl), numeric(1))
    crown <- max(abs(mean_fit - mean_exact) / mean_exact)
    cat(sprintf("%-24s %5d %10.3e %10.3e %10.3e %10.3e %10.3e %11.3e\n",
                nm, length(kx), ev[["max"]], ev[["rms"]], es[["max"]],
                es[["rms"]], min(fv$value), crown))
  }

  cat(sprintf("\n%-24s %5s %12s %12s\n", "dropped position channel", "knots",
              "max/range", "rms/range"))
  for (nm in names(cands)) {
    place <- cands[[nm]]
    d <- tryCatch(dropped_channel(patch, place, hmax, z, space),
                  error = function(e) c(max = NA, rms = NA))
    cat(sprintf("%-24s %5d %12.3e %12.3e\n", nm, length(place(hmax)),
                d[["max"]], d[["rms"]]))
  }
  invisible(NULL)
}

study(ladder_patch_two_by_two(), "open stand, 2 species x 2 cohorts")
study(shaded_patch(),            "shaded stand, 8 cohorts, L(0)=0.055")
