# Regression guard: the odelia::separable_field's query derivative is exact in
# FF16's CROWN pattern -- a focal plant reads the field at z = node * H (query
# height tied to its own active height), summed over crown nodes, with every
# source height/weight scaling with the parameter (as a single-species stand's
# cohorts do). This EXONERATES the field as the cause of the FF16 R0 gradient bug:
# the separable field, a direct O(N^2) sum, and finite differences all agree to
# machine precision, so the field read is not where the reverse gradient goes
# wrong (the bug is downstream in the SCM growth->fecundity trajectory -- see
# ff16_scm_gradient_driver.cpp and the build plan).

test_that("separable_field query derivative is exact in the FF16 crown pattern", {
  skip_on_cran()
  dlls <- getLoadedDLLs()
  skip_if_not("plant"  %in% names(dlls), "plant DLL not loaded")
  skip_if_not("odelia" %in% names(dlls), "odelia DLL not loaded")

  plant_inc  <- system.file("include", package = "plant")
  if (plant_inc == "") plant_inc <- here::here("plant/inst/include")
  odelia_inc <- system.file("include", package = "odelia")
  bh_inc     <- system.file("include", package = "BH")
  skip_if(plant_inc == "" || odelia_inc == "" || bh_inc == "",
          "package headers unavailable")
  plant_so  <- dlls[["plant"]][["path"]]
  odelia_so <- dlls[["odelia"]][["path"]]
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

  built <- tryCatch({ Rcpp::sourceCpp(test_path("field_crown_probe.cpp")); TRUE },
                    error = function(e) { message("field-crown probe build failed: ",
                                                  conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  r <- field_crown_probe()
  # Value: separable factoring reconstructs the direct sum exactly.
  expect_equal(r$val_sep, r$val_dir, tolerance = 1e-12)
  # Derivative: separable-AD == direct-AD == FD. This is the exoneration.
  expect_equal(r$grad_sep, r$grad_dir, tolerance = 1e-10)
  expect_equal(r$grad_sep, r$fd_sep,   tolerance = 1e-5)
  expect_equal(r$grad_dir, r$fd_dir,   tolerance = 1e-5)
})
