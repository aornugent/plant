# Stage E (light channel, ad-implementation.md §4.2 M3 / §7.3): the resident
# self-shading light -> leaf coupling. The canopy openness the leaf sees is
# active on a resident pass (the light field is built from the stand's own
# competition), so the leaf supplied_derivative seam injects d(profit)/d(light).
# This also carries k_I through the leaf (radiation = k_I * openness * PPFD, now
# recomputed inside leaf_profit_frozen).
#
# This test seeds the resident light active, computes the growth rate through the
# leaf at one compute_rates, and FD-checks d(growth)/d(light) for both single-
# solve shading models (mean-light -- TF24's default -- and crown-centre).
#
# Compiles a driver with sourceCpp; skips where that is unavailable.

test_that("TF24 resident light->leaf gradient matches finite difference", {
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
    Rcpp::sourceCpp(test_path("tf24_light_coupling_driver.cpp"))
    TRUE
  }, error = function(e) { message("TF24 light-coupling sourceCpp build failed: ",
                                   conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  for (sm in c("", "crown-centre")) {          # "" == mean-light (default)
    for (lt in c(0.4, 0.6, 0.9)) {
      r <- tf24_light_channel_check(shading_model = sm, light = lt,
                                    height = 1.0, delta = 1e-6)
      expect_equal(r$ad_grad, r$fd_grad, tolerance = 1e-4,
                   info = sprintf("d(growth)/d(light) [%s] at light=%.1f: ad=%.8g fd=%.8g",
                                  ifelse(sm == "", "mean-light", sm), lt,
                                  r$ad_grad, r$fd_grad))
    }
  }
})
