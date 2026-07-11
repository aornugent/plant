# The FF16 resident self-shading census gradient through the compiled plant.so
# path. stand_gradient recomputes the canopy live on the recorded knot positions
# with the active cohorts (L3 cache empty), so a trait re-shades the stand through
# area_leaf and the self-shading cross term flows -- the coupling #540 mocked with a
# hand-built Beer's-law spline now running as the real environment recompute. The
# recorded baseline is the tight regression oracle; the self-shading term is shown
# present (the resident gradient differs from the frozen invasion one) and bracketed
# by a full-run FD. AD needs the installed DLL, so this skips under load_all.

test_that("FF16 resident self-shading gradient runs on the compiled path", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  scm <- gradient_fixture_scm()
  traits <- gradient_fixture_traits
  metrics <- gradient_fixture_metrics$resident

  res <- stand_gradient(scm, metrics, traits)
  expect_equal(dim(res$gradient), c(length(metrics), length(traits)))
  expect_true(all(is.finite(res$gradient)))

  base <- read_gradient_baseline()
  expect_matches_gradient_baseline(res$gradient, base$jacobians$resident$gradient,
                                   base$fingerprint, "resident census gradient")
})

test_that("the resident self-shading term is present (differs from the frozen invasion gradient)", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  scm <- gradient_fixture_scm()
  traits <- gradient_fixture_traits
  metrics <- gradient_fixture_metrics$resident

  res <- stand_gradient(scm, metrics, traits)$gradient
  inv <- invasion_gradient(scm, metrics, traits)$gradient[metrics, traits]

  # Recomputing the canopy adds the self-shading cross term the frozen invasion
  # gradient drops, so the two Jacobians differ materially on the dominant column.
  expect_gt(max(abs(res[, "lma"] - inv[, "lma"])), 1)
})

test_that("the resident gradient is bracketed by a full self-shading FD", {
  skip_if(is_pkgload_dll_plant(), "AD oracle recorded against the installed plant.so")

  pr <- gradient_fixture_parameters()
  e <- Environment("FF16")
  ctrl <- Control(); ctrl$save_RK45_cache <- TRUE
  scm <- run_scm(pr, e, ctrl)
  traits <- gradient_fixture_traits
  metrics <- gradient_fixture_metrics$resident

  ad <- stand_gradient(scm, metrics, traits)$gradient

  # Central FD of the FULL self-shading run: perturb the trait and re-run the whole
  # resident. Unlike the frozen-canopy invasion FD, this moves the adaptive knot
  # placement with the trait, so it also differences the node-placement derivative
  # the fixed-knot replay holds constant by design (spec S1). AD is the exact
  # derivative of the recorded replay; the gap to this FD is that accepted
  # record->replay residual, so the check is a same-sign order-of-magnitude bracket,
  # not the 1e-4 the clean invasion axis meets. The recorded baseline is the tight
  # regression oracle.
  fd_val <- function(param, delta) {
    p <- pr
    p$strategies[[1]]$pars[[param]] <- p$strategies[[1]]$pars[[param]] + delta
    stand_gradient(run_scm(p, e, ctrl), metrics, param)$value
  }
  central <- function(param, x0, h) (fd_val(param, h * x0) - fd_val(param, -h * x0)) / (2 * h * x0)

  for (tr in traits) {
    x0 <- pr$strategies[[1]]$pars[[tr]]
    fd <- central(tr, x0, 1e-4)
    for (mt in metrics) {
      a <- unname(ad[mt, tr]); f <- unname(fd[[mt]])
      expect_gt(a * f, 0)                                  # same sign
      expect_lt(abs(a - f) / max(abs(f), 1e-30), 0.15)     # within the replay residual
    }
  }
})
