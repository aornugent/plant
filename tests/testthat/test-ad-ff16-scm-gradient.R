# FF16 SCM R0 (offspring) gradient over the whole FF16 method-of-characteristics
# solve. The active pass replays the RESOLVED schedule (L0 node_schedule_times + L1
# ode_times from run_scm(refine_schedule=TRUE)), i.e. exactly what
# run_scm(use_ode_times=TRUE) replays.
#
# STATUS (2026-07-20): FF16 now transports on the geometric mass chart (lambda =
# log density + log spacing, evolving by -loss with the compression cancelled
# identically; K93's machinery, gated by the geometric_transport marker + the
# node_geometric_compression Control flag). Three settled facts.
#  (1) NO schedule sensitivity. The frozen replay on the resolved schedule matches
#      the fully-adaptive model in double. The correct replay is the resolved L0+L1.
#  (2) The chart made the pure-loss channel EXACT: d(offspring)/d(k_l) reverse AD
#      == resolved-schedule FD to ~1e-9 (was broken before the chart). k_l enters
#      as loss only, so it flows straight through the -loss transport with no
#      compression residual. lma tightened to ~0.36% (dominant magnitude).
#  (3) THE ADJOINT GAP IS CLOSED (superseded 2026-07-27). This header used to record a
#      residual on the growth+self-shading traits -- a_l1 reverse AD 0.0629 vs FD
#      plateau 0.1007, a 38% gap -- and that is NO LONGER TRUE. The a1-a4 severance
#      fixes (templated CanopyShape, smooth_positive on the growth/fecundity clamp,
#      the IFT-lifted birth size) closed it. Measured on this file's own driver:
#        lma  -56.0698330 vs FD -56.0699810   rel 2.64e-06
#        a_l1   0.1007215 vs FD   0.1007221   rel 6.06e-06
#        k_l   -2.6635231 vs FD  -2.6635231   rel 1.31e-08
#      So all three channels agree to the FD's own noise floor. NOTE the assertions
#      below still gate lma and a_l1 at 1e-2, which is ~4 orders looser than the truth
#      -- a 100x accuracy regression would pass green. Tightening them is deliberate
#      future work, not an oversight to fix blindly: verify the FD reference first
#      (HANDOFF Part 1), since these are the channels whose FD is delta-sensitive.
#      Recorded in docs/v3-requirements.md 7.1 and 6/B5.
#
# life = 40: FF16 has reproduced (offspring ~0.6) and the reverse tape fits (life 50
# exceeds memory; checkpointing deferred).

test_that("FF16 SCM offspring (R0) gradient: value exact; loss channel exact; self-shading gap known", {
  skip_on_cran()
  dlls <- getLoadedDLLs()
  skip_if_not("plant"  %in% names(dlls), "plant DLL not loaded")
  skip_if_not("odelia" %in% names(dlls), "odelia DLL not loaded")

  plant_so  <- dlls[["plant"]][["path"]]
  odelia_so <- dlls[["odelia"]][["path"]]
  plant_inc  <- system.file("include", package = "plant")
  if (plant_inc == "") plant_inc <- here::here("plant/inst/include")
  odelia_inc <- system.file("include", package = "odelia")
  bh_inc     <- system.file("include", package = "BH")
  skip_if(plant_inc == "" || odelia_inc == "" || bh_inc == "",
          "package headers unavailable")
  skip_if(!file.exists(plant_so) || !file.exists(odelia_so),
          "compiled plant/odelia libraries unavailable")

  old <- Sys.getenv(c("PKG_CXXFLAGS", "PKG_LIBS"), unset = NA)
  Sys.setenv(PKG_CXXFLAGS = paste0("-isystem", plant_inc, " -I", odelia_inc,
                                   " -I", bh_inc))
  Sys.setenv(PKG_LIBS = paste(shQuote(odelia_so), shQuote(plant_so)))
  on.exit({
    for (k in names(old)) {
      if (is.na(old[[k]])) Sys.unsetenv(k)
      else do.call(Sys.setenv, setNames(list(old[[k]]), k))
    }
  }, add = TRUE)

  built <- tryCatch({ Rcpp::sourceCpp(test_path("ff16_scm_gradient_driver.cpp")); TRUE },
                    error = function(e) { message("FF16 gradient build failed: ",
                                                  conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  lma0 <- 0.1978791; life <- 40
  p0 <- scm_base_parameters("FF16"); p0$max_patch_lifetime <- life
  pr <- add_strategies(p0, trait_matrix(lma0, "lma"), birth_rate = 20)
  base <- run_scm(pr, refine_schedule = TRUE)          # resolve L0 + L1
  nst <- base$parameters$node_schedule_times
  ot  <- base$parameters$ode_times

  r <- ff16_scm_offspring_gradient(nst, ot, lma = lma0, max_patch_lifetime = life)

  # Value: the active replay reproduces the double offspring on the resolved schedule.
  expect_equal(r$value, r$offspring_double, tolerance = 1e-6)
  expect_true(r$value > 0)

  # Self-consistency smoke test (necessary, NOT sufficient): reverse == forward.
  expect_equal(r$jvp, r$dot, tolerance = 1e-6)

  # Target order is {lma (0), a_l1 (6), k_l (16)} = FF16_R0_TARGET_IDX.
  # Correctness gate, per channel:
  #  - k_l (loss only): the mass chart transports loss exactly, so reverse AD
  #    matches the resolved-schedule FD to numerical tolerance. This is the chart's
  #    load-bearing win; it turns red if a regression re-severs the loss adjoint.
  expect_equal(r$grad[[3]], r$fd_grad[[3]], tolerance = 1e-6)
  #  - lma (growth + self-shading): within ~0.5% of FD; a small self-shading
  #    residual remains (see below), but lma's magnitude keeps it tight.
  expect_equal(r$grad[[1]], r$fd_grad[[1]], tolerance = 1e-2)
  #  - a_l1 (growth + self-shading): matches FD once the growth/fecundity clamp is
  #    smoothed (util::smooth_positive) and the birth-size channel is lifted
  #    (Patch::reset re-prepares the strategy from the seeded parameters, and
  #    area_leaf_0 is derived from the lifted initial_height_). The earlier gap was
  #    the hard-clamp kink plus the severed birth size, not a self-shading residual.
  expect_equal(r$grad[[2]], r$fd_grad[[2]], tolerance = 1e-2)
})
