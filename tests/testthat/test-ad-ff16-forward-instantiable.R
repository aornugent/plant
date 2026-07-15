# The FF16 rate path is forward-mode instantiable: its rebind alias lets generic
# code lift the whole strategy/environment to a nested forward-over-reverse type,
# which is what the Species-level geometric compression gates on (census
# gradients). This test compiles a driver that instantiates Species<FF16_active>
# at both the active reverse scalar and the nested forward type, then runs
# compute_rates through the geometric and stencil branches. It is primarily a
# COMPILE guard -- the assimilation quadrature (QK), the light spline, and the
# birth-height lift each narrow to double in places that must strip every AD
# layer, or the build fails.
#
# Like the Gate 0 driver, it links the compiled plant + odelia shared libraries
# and skips (rather than fails) where sourceCpp cannot build in the session.

test_that("FF16 rate path instantiates at active + nested forward-over-reverse", {
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
    Rcpp::sourceCpp(test_path("ff16_forward_instantiable_driver.cpp"))
    TRUE
  }, error = function(e) { message("FF16 forward-instantiable sourceCpp build ",
                                   "failed: ", conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  r <- ff16_forward_instantiable_smoke()

  # No parameter is seeded, so every path shares the same primal: the active
  # and double runs, and the geometric and stencil branches (which differ only
  # in the log-density transport rate, not in height_max), must all agree.
  expect_equal(r$height_max_active_stencil,   r$height_max_double_stencil)
  expect_equal(r$height_max_active_geometric, r$height_max_double_geometric)
  expect_equal(r$height_max_double_stencil,   r$height_max_double_geometric)
  expect_true(is.finite(r$height_max_double_stencil))
})
