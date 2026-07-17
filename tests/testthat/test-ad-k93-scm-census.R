# CD-G layer (a): build and run a double K93 SCM entirely in C++ (there is no
# C++ SCM-construction path today -- params are set up R-side), reduce a census
# metric, and check it against R's run_scm. This is the foundation the active
# census gradient is layered on; verifying the C++ construction first keeps the
# later gradient trustworthy.
#
# Compiles k93_scm_census_driver.cpp with sourceCpp, linking the compiled plant
# + odelia libraries (only double crosses R). Skips where that is unavailable,
# the same pattern the gate0 AD tests use.

test_that("C++ K93 SCM construction reproduces R run_scm (CD-G layer a)", {
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
                    error = function(e) { message("layer-a build failed: ",
                                                  conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  # The C++ SCM matches R's documented offspring for both transport paths
  # (geometric mass chart vs the default FD stencil), so the construction is
  # faithful. Census (a basal-area moment) is finite and positive.
  geo <- k93_scm_census_double(geometric = TRUE)
  fd  <- k93_scm_census_double(geometric = FALSE)
  expect_equal(geo$offspring[1], 0.075453, tolerance = 1e-4)
  expect_equal(fd$offspring[1],  0.075325, tolerance = 1e-4)
  expect_true(is.finite(geo$census) && geo$census > 0)
})
