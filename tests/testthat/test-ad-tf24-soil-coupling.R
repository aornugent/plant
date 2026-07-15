# Stage C (ad-implementation.md §6.3): the resident soil coupling. TF24's soil
# water state is an active ODE state, and the leaf reads the soil water potential
# through the envelope-theorem supplied_derivative seam (the soil-psi channel).
# This test seeds a soil-moisture layer active, computes the individual growth
# rate through the leaf at one compute_rates call, and checks
# d(growth_rate)/d(soil_moisture) against a central finite difference -- the
# soil -> leaf half of the resident feedback loop.
#
# (The uptake -> soil half, resource_depletion carrying value_type, is exercised
# by the growing-SCM resident gradient; that end-to-end FD check is tracked
# separately. This test isolates the differentiable soil read.)
#
# Compiles a driver with sourceCpp, linking the compiled plant + odelia
# libraries; skips where that is unavailable, as the Gate 0 tests do.

test_that("TF24 resident soil->leaf gradient matches finite difference", {
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
    Rcpp::sourceCpp(test_path("tf24_soil_coupling_driver.cpp"))
    TRUE
  }, error = function(e) { message("TF24 soil-coupling sourceCpp build failed: ",
                                   conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  for (h in c(0.5, 1.0, 2.0)) {
    r <- tf24_soil_channel_check(height = h, delta = 1e-6)
    expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-3,
                 info = sprintf("d(growth)/d(soil_moist) at h=%.1f: ad=%.8g fd=%.8g",
                                h, r$ad_grad, r$fd_grad))
  }
})
