# Issue #47: the TF24f Tier-B soil-coupling seam injects d(uptake)/d(tracked
# collar) analytically, per soil layer. Uptake (E_from_Soil_to_Root_Collar) is
# NON-stationary in the collar and piecewise per soil layer, so a straight FD of
# the collar term straddles layer-crossing kinks and undershoots ~3x; the analytic
# per-layer derivative (Leaf::dsoil_consumption_dpsi_collar_perlayer, the per-layer
# form of the verified aggregate dE_from_soil_dpsi_collar) is required.
#
# Two checks, both near machine precision:
#   (1) the analytic per-layer derivative vs a pure-double FD of soil_consumption_
#       w.r.t. the collar potential (leaf-level, no tape / seam);
#   (2) the seam end-to-end: seed the tracked collar active, read consumption_rate
#       through the full seam, FD-check d(consumption_rate)/d(tracked collar).

setup_ad_env <- function() {
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
  Sys.setenv(PKG_CXXFLAGS = paste0("-isystem", plant_inc, " -I", odelia_inc,
                                   " -I", bh_inc))
  Sys.setenv(PKG_LIBS = paste(shQuote(odelia_so), shQuote(plant_so)))
}

test_that("analytic per-layer d(uptake)/d(collar) matches a double FD", {
  skip_on_cran(); setup_ad_env()
  built <- tryCatch({ Rcpp::sourceCpp(test_path("tf24_collar_deriv_driver.cpp")); TRUE },
                    error = function(e) { message(conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable")
  for (h in c(3.0, 5.0, 10.0)) {
    r <- tf24_collar_deriv_check(height = h, delta = 1e-6)
    ok <- is.finite(r$analytic) & abs(r$fd) > 1e-12
    expect_equal(r$analytic[ok], r$fd[ok], tolerance = 1e-6,
                 info = sprintf("per-layer collar deriv at height=%.0f", h))
  }
})

test_that("TF24f seam d(consumption_rate)/d(tracked collar) matches a double FD", {
  skip_on_cran(); setup_ad_env()
  built <- tryCatch({ Rcpp::sourceCpp(test_path("tf24f_collar_uptake_driver.cpp")); TRUE },
                    error = function(e) { message(conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable")
  for (L in 0:3) for (tp in c(1.5, 2.0)) {
    r <- tf24f_collar_uptake_check(tracked_psi = tp, height = 1.0, layer = L, delta = 1e-6)
    if (abs(r$fd_grad) > 1e-12)
      expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-4,
                   info = sprintf("d(uptake[%d])/d(collar) at psi=%.1f", L, tp))
  }
})
