# Gate 0 for TF24 (ad-implementation.md §15, §7.3): the TF24 single-plant
# IndividualRunner is the first active reverse-mode AD target for TF24. Unlike
# FF16/K93, TF24's growth runs an embedded double Leaf hydraulic optimiser; its
# parameter sensitivity reaches the tape through the envelope-theorem
# supplied_derivative seam (d(profit)/d(theta) by central FD of the leaf profit
# at the FROZEN optimum collar psi). This test seeds one low-level TF24 parameter,
# replays a fixed schedule at an active scalar, and checks the AD gradient of
# final plant height against a central finite difference.
#
# It compiles a small driver (gate0_tf24_driver.cpp) with sourceCpp, linking the
# compiled plant + odelia shared libraries. Where that is unavailable (e.g. a
# stripped CI image) the test skips rather than fails, as the FF16 Gate 0 does.

test_that("TF24 IndividualRunner AD gradient matches finite difference (Gate 0)", {
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
    Rcpp::sourceCpp(test_path("gate0_tf24_driver.cpp"))
    TRUE
  }, error = function(e) { message("TF24 Gate 0 sourceCpp build failed: ",
                                   conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  # Leaf-hydraulic parameters: they act ONLY through the envelope-theorem leaf
  # seam (no birth-height channel), so the FD oracle is clean and the injected
  # partial matches it to the leaf FD floor.
  for (p in c("vcmax_25", "jmax_25", "K_s", "k_I")) {
    r <- tf24_gate0_fd_check(p, t_end = 5.0, delta = 1e-5)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-3,
                 info = sprintf("d(height)/d(%s): ad=%.8g fd=%.8g", p,
                                r$ad_grad, r$fd_grad))
  }

  # Allometry / mass-cascade parameters that also move the birth height h*(theta)
  # through mass_live_given_height. h* is a double root-find whose parameter
  # derivative is lifted onto the tape by the implicit function theorem
  # (lift_birth_height); the AD value is exact and the looser tolerance reflects
  # the FD ORACLE, which re-solves h* per perturbation (a_l1 is birth-dominant).
  for (p in c("lma", "rho", "a_l1")) {
    r <- tf24_gate0_fd_check(p, t_end = 5.0, delta = 1e-4)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1.5e-2,
                 info = sprintf("d(height)/d(%s): ad=%.8g fd=%.8g", p,
                                r$ad_grad, r$fd_grad))
  }
})
