# Gate 0 (ad-implementation.md §15): the FF16 single-plant IndividualRunner is
# the first active reverse-mode AD target. This test seeds one low-level FF16
# parameter, replays a fixed schedule at an active scalar, and checks the AD
# gradient of final plant height against a central finite difference.
#
# It compiles a small driver (gate0_ff16_driver.cpp) with sourceCpp, linking the
# compiled plant + odelia shared libraries for the XAD Tape / QK runtime symbols
# (only double crosses R, so the active types never leave C++). The build needs
# the package headers and both DLLs loaded; where that is unavailable (e.g. a
# stripped CI image), the test skips rather than fails -- the same pattern
# odelia's sourceCpp AD tests use.

test_that("FF16 IndividualRunner AD gradient matches finite difference (Gate 0)", {
  skip_on_cran()
  dlls <- getLoadedDLLs()
  skip_if_not("plant"  %in% names(dlls), "plant DLL not loaded")
  skip_if_not("odelia" %in% names(dlls), "odelia DLL not loaded")

  plant_so   <- dlls[["plant"]][["path"]]
  odelia_so  <- dlls[["odelia"]][["path"]]
  plant_inc  <- system.file("include", package = "plant")
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

  built <- tryCatch({
    Rcpp::sourceCpp(test_path("gate0_ff16_driver.cpp"))
    TRUE
  }, error = function(e) { message("Gate 0 sourceCpp build failed: ",
                                   conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  # Parameters whose derivative is fully on the taped rate path (directly, or via
  # a prepare_strategy-cached S quantity that reset() recomputes, e.g. eta ->
  # eta_c). These must match to well within FD accuracy.
  for (p in c("a_p1", "a_p2", "k_l", "eta")) {
    r <- ff16_gate0_fd_check(p, t_end = 5.0, delta = 1e-5)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-4,
                 info = sprintf("d(height)/d(%s): ad=%.8g fd=%.8g", p,
                                r$ad_grad, r$fd_grad))
  }

  # lma / a_l1 additionally act through the birth height (height_0), which is a
  # double root-find until the height_seed supplied_derivative seam (§7.2) lands.
  # AD captures the rate-path part; the birth-height part is the known Stage-B
  # residual, so they are NOT expected to match yet -- asserting the machinery
  # runs for them (finite, correct sign) without demanding the deferred term.
  r_lma <- ff16_gate0_fd_check("lma", t_end = 5.0, delta = 1e-5)
  expect_true(is.finite(r_lma$ad_grad))
  expect_lt(r_lma$ad_grad, 0)  # more leaf mass per area slows height growth
})
