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
  expect_equal(r$value, r$value_double, tolerance = 1e-10)

  # Self-consistency smoke test (necessary, NOT sufficient): reverse vs forward.
  # This alone once gave false confidence -- both legs shared the frozen-query
  # bias -- so it is no longer the correctness gate.
  expect_equal(r$jvp, r$dot, tolerance = 1e-6)

  # The correctness gate: the reverse gradient matches a pinned-schedule central
  # FD of the model (catches the shared-bias class the oracle cannot). The exact
  # separable_field makes the self-shading feedback flow, so these now agree.
  expect_equal(r$grad, r$fd_grad, tolerance = 1e-3)

  # Structure: b_0, b_1 (growth) drive the basal-area census; d_0 (recruitment)
  # has no path to it -- a structural zero the oracle cannot check.
  expect_true(abs(r$grad[1]) > 1)      # d(census)/d(b_0)
  expect_true(abs(r$grad[2]) > 1)      # d(census)/d(b_1)
  expect_equal(r$grad[3], 0)           # d(census)/d(d_0): structural zero
})

test_that("K93 SCM offspring (R0) gradient passes the JVP=VJP oracle (CD-G layer b)", {
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

  # The offspring (R0) functional is differentiated through the value_type-carrying
  # reproduction chain (node fecundity -> species -> patch -> trapezium). The active
  # value reproduces the double offspring production exactly.
  o <- k93_scm_offspring_gradient()
  expect_equal(o$value, o$value_double, tolerance = 1e-10)
  expect_equal(o$value, 0.0754715, tolerance = 1e-4)  # documented K93 offspring (mass chart)

  # Self-consistency smoke test (necessary, not sufficient).
  expect_equal(o$jvp, o$dot, tolerance = 1e-6)
  # The correctness gate: reverse gradient vs the pinned-schedule model FD.
  expect_equal(o$grad, o$fd_grad, tolerance = 1e-3)

  # Targets {b_0, d_0, d_1}. Unlike the census, all three feed the fitness
  # integral: d_0, d_1 are the recruitment/suppression traits (the dominant
  # contributors) and b_0 is growth. This is the complement of the census's
  # structural zero -- the reproduction chain carries derivatives end to end.
  expect_true(abs(o$grad[2]) > 1)      # d(R0)/d(d_0): recruitment, dominant
  expect_true(o$grad[2] > 0)           # more recruitment -> more offspring
  expect_true(abs(o$grad[3]) > 0)      # d(R0)/d(d_1): suppression, nonzero
})
