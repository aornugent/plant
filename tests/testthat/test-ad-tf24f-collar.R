# Stage D (issue #43): TF24f's leaf runs at a TRACKED collar-psi ODE state, not
# the optimum, so the envelope theorem does not apply -- d(profit)/d(psi) != 0.
# The leaf supplied_derivative seam gains a collar-psi channel that injects
# d(profit)/d(psi) (the leaf's analytic gradient, the same one the acclimation
# rate uses) onto the active tracked state, so the acclimation channel
# d(profit)/d(psi) * d(psi_tracked)/d(theta) is on the tape.
#
# This test seeds the tracked collar-psi state active (off the optimum), computes
# the growth rate through the leaf at one compute_rates call, and checks
# d(growth)/d(tracked_psi) against a central finite difference. The injected
# partial is analytic, so agreement is near machine precision.
#
# Compiles a driver with sourceCpp; skips where that is unavailable.

test_that("TF24f tracked collar-psi gradient matches finite difference", {
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
    Rcpp::sourceCpp(test_path("tf24f_collar_driver.cpp"))
    TRUE
  }, error = function(e) { message("TF24f collar sourceCpp build failed: ",
                                   conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  for (psi in c(0.8, 1.5, 2.5)) {
    r <- tf24f_collar_check(tracked_psi = psi, height = 1.0, delta = 1e-6)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-4,
                 info = sprintf("d(growth)/d(tracked_psi) at psi=%.1f: ad=%.8g fd=%.8g",
                                psi, r$ad_grad, r$fd_grad))
  }
})
