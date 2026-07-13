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

  # Parameters that act ONLY on the taped rate path -- no inner solve on the
  # metric's graph, so the FD oracle is itself exact and AD matches it to machine
  # precision (~1e-11).
  for (p in c("a_p1", "a_p2", "k_l")) {
    r <- ff16_gate0_fd_check(p, t_end = 5.0, delta = 1e-4)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-4,
                 info = sprintf("d(height)/d(%s): ad=%.8g fd=%.8g", p,
                                r$ad_grad, r$fd_grad))
  }

  # Parameters that also act through the birth height h*(theta): a different
  # parameter gives a different-sized seedling (lma and a_l1 through the allometry;
  # eta through eta_c in the sapwood/bark mass). h* is a double root-find whose
  # parameter derivative is lifted onto the tape by the implicit function theorem
  # (FF16_Strategy_::lift_birth_height). The AD value is EXACT and invariant to the
  # FD step; the looser tolerance reflects the FD ORACLE, which re-solves h* per
  # perturbation, so a small delta is dominated by its ~1e-8 root noise (see the
  # FD-verification note in ad-implementation.md §15). delta=1e-3 sits in the
  # oracle's clean band for all three; a delta-sweep confirms FD -> AD there.
  for (p in c("lma", "a_l1", "eta")) {
    r <- ff16_gate0_fd_check(p, t_end = 5.0, delta = 1e-3)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-3,
                 info = sprintf("d(height)/d(%s): ad=%.8g fd=%.8g", p,
                                r$ad_grad, r$fd_grad))
  }
})
