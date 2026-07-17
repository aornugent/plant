# FF16 SCM R0 (offspring) gradient: the exact reverse-mode trait gradient of
# lifetime offspring production over the whole FF16 method-of-characteristics
# solve, verified by the JVP=VJP dot-product oracle.
#
# This is the second strategy on the v2 SCM-gradient stack after K93, and it
# reuses the generic machinery with no strategy-specific gradient code: the SCM
# duck-types the odelia gradient-driver contract, the active pass replays the L1
# ode-time schedule recorded on a double run (resident L2, no set_mutant), and
# offspring flows through the value_type reproduction chain. The only FF16-side
# work was making its competition/environment path active-instantiable (the ratio
# collapse in compute_competition_by_ratio and building the light field at the
# active scalar) -- both bit-identical on the double path.
#
# Uses a shortened max_patch_lifetime (50): the full lifetime (105.32) produces a
# reverse tape that exceeds memory for FF16's heavy rate path (leaf photosynthesis
# quadrature + allocation + birth-height IFT per cohort per step) -- full-lifetime
# runs need tape checkpointing at the node-introduction boundary (deferred). At 50,
# FF16 has matured and reproduced (offspring ~3.4, well clear of zero), so the
# oracle is a real test, not a trivial 0 = 0.

test_that("FF16 SCM offspring (R0) gradient passes the JVP=VJP oracle", {
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

  r <- ff16_scm_offspring_gradient(max_patch_lifetime = 50)

  # The active replay reproduces the double offspring closely (~1e-9). It is not
  # bit-exact: FF16 rebuilds its light-availability spline at the active scalar
  # (the L2 recompute), whose adaptive knot values differ from the double field at
  # the ~1e-9 level. Physically negligible, and the oracle below independently
  # confirms the AD is self-consistent regardless.
  expect_equal(r$value, r$offspring_double, tolerance = 1e-6)
  # A real (non-trivial) fitness -- FF16 has matured and reproduced.
  expect_true(r$value > 1)

  # The FD-free correctness gate.
  expect_equal(r$jvp, r$dot, tolerance = 1e-6)

  # Targets {lma, a_l1, k_l}. lma (leaf mass per area) is the classic FF16
  # cost trait: heavier leaves lower lifetime fitness, so d(R0)/d(lma) < 0.
  expect_true(r$grad[1] < 0)
})
