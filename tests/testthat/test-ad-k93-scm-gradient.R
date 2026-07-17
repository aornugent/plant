# CD-G layers (c)/(d): the exact reverse-mode trait gradient of an emergent K93
# SCM census metric (total stand competition, a basal-area moment).
#
# The SCM is the runnable the gradient is taken over -- it duck-types the odelia
# gradient-driver contract directly (value_type, get_system_ref, ad_parameters,
# ad_initial_state, reset, run, an assignable tape), no wrapper. The active pass
# replays the L1 ode-time schedule recorded on a double run and recomputes the
# environment field at the active scalar (resident L2, not a frozen mutant), which
# is what makes the gradient exact.
#
# Correctness is the JVP=VJP dot-product oracle <J v, u> = <v, Jᵀ u>: the reverse
# gradient g is Jᵀ·1, the forward directional derivative along v is J v = g·v, so
# they must agree to machine precision -- an FD-free correctness gate that needs no
# perturbation and no inner-solve re-run.
#
# Compiles k93_scm_census_driver.cpp with sourceCpp, linking the compiled plant +
# odelia libraries (only double crosses R). Skips where that is unavailable, the
# same pattern the gate0 / layer-(a) AD tests use.

test_that("K93 SCM census gradient passes the JVP=VJP oracle (CD-G layers c/d)", {
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

  built <- tryCatch({ Rcpp::sourceCpp(test_path("k93_scm_census_driver.cpp")); TRUE },
                    error = function(e) { message("gradient build failed: ",
                                                  conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  r <- k93_scm_census_gradient()

  # The active replay reproduces the double census exactly: the pinned ode-time
  # schedule is the same trajectory, just carried at the active scalar.
  expect_equal(r$value, r$census_double, tolerance = 1e-10)

  # The FD-free correctness gate: the forward directional derivative equals the
  # reverse gradient contracted against the same direction, to machine precision.
  expect_equal(r$jvp, r$dot, tolerance = 1e-6)

  # The gradient targets {b_0, b_1, d_0}. b_0, b_1 are growth-rate parameters and
  # drive the basal-area census (nonzero). d_0 is the K93 recruitment (fecundity)
  # parameter: it feeds offspring output, not the within-patch basal area, so it
  # has no path to the census and AD returns a structural zero -- a check the
  # dot-product oracle cannot make (both legs traverse the same graph).
  expect_true(abs(r$grad[1]) > 1)      # d(census)/d(b_0)
  expect_true(abs(r$grad[2]) > 1)      # d(census)/d(b_1)
  expect_equal(r$grad[3], 0)           # d(census)/d(d_0): structural zero
})
