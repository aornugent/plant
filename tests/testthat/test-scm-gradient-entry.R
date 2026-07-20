# The run-shaped gradient entry (plant/scm_gradient.h): scm_gradient/scm_jacobian.
#
# The entry OWNS record->replay. A caller passes only trait values + target field
# indices + a functional; the entry runs the adaptive double refine itself, crosses
# double->active via the odelia System rebind_from contract, and replays the resolved
# schedule. There is no schedule argument, so a caller cannot express a wrong replay
# grid (the -255/+442 step_history footgun the bespoke drivers had is unrepresentable).
#
# Gate: the entry must (1) reproduce the double value exactly on its self-resolved
# schedule, and (2) return the SAME certified gradient as the bespoke driver -- which
# is handed R's run_scm(refine_schedule) schedule -- proving the entry's internal
# refine reproduces run_scm's, and match that driver's pinned-schedule model FD to the
# certified per-trait tolerances.

setup_link <- function() {
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
  Sys.setenv(PKG_CXXFLAGS = paste0("-isystem", plant_inc, " -I", odelia_inc,
                                   " -I", bh_inc))
  Sys.setenv(PKG_LIBS = paste(shQuote(odelia_so), shQuote(plant_so)))
}

test_that("FF16 entry gradient: self-refines the schedule, reproduces the certified driver", {
  skip_on_cran()
  old <- Sys.getenv(c("PKG_CXXFLAGS", "PKG_LIBS"), unset = NA)
  setup_link()
  on.exit({
    for (k in names(old)) {
      if (is.na(old[[k]])) Sys.unsetenv(k)
      else do.call(Sys.setenv, setNames(list(old[[k]]), k))
    }
  }, add = TRUE)

  built <- tryCatch({
    Rcpp::sourceCpp(test_path("scm_gradient_driver.cpp"))
    Rcpp::sourceCpp(test_path("ff16_scm_gradient_driver.cpp"))
    TRUE
  }, error = function(e) { message("build failed: ", conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  lma0 <- 0.1978791; life <- 40

  # The entry: no schedule in, it refines internally.
  ent <- ff16_entry_gradient(c(0L, 6L, 16L), lma = lma0, birth_rate = 20,
                             max_patch_lifetime = life, metric = 0L)

  # The certified bespoke driver: handed R's run_scm(refine_schedule) schedule.
  p0 <- scm_base_parameters("FF16"); p0$max_patch_lifetime <- life
  pr <- add_strategies(p0, trait_matrix(lma0, "lma"), birth_rate = 20)
  base <- run_scm(pr, refine_schedule = TRUE)
  bespoke <- ff16_scm_offspring_gradient(base$parameters$node_schedule_times,
                                         base$parameters$ode_times,
                                         lma = lma0, max_patch_lifetime = life)

  # (1) The active replay reproduces the double value on the entry's own schedule.
  expect_equal(ent$value, ent$value_double, tolerance = 1e-10)

  # (2) The entry's self-refined schedule reproduces run_scm's: same value and same
  #     reverse gradient as the bespoke driver, to numerical tolerance.
  expect_equal(ent$value, bespoke$value, tolerance = 1e-6)
  expect_equal(ent$grad, bespoke$grad, tolerance = 1e-6)

  # (3) The certified per-trait model FD (bespoke driver's pinned-schedule FD):
  #     k_l (loss only) exact; lma / a_l1 (growth + self-shading) within ~1%.
  expect_equal(ent$grad[[3]], bespoke$fd_grad[[3]], tolerance = 1e-5)
  expect_equal(ent$grad[[1]], bespoke$fd_grad[[1]], tolerance = 1e-2)
  expect_equal(ent$grad[[2]], bespoke$fd_grad[[2]], tolerance = 1e-2)
})

test_that("K93 entry gradient: self-refines the schedule, matches a reoptimising FD", {
  skip_on_cran()
  old <- Sys.getenv(c("PKG_CXXFLAGS", "PKG_LIBS"), unset = NA)
  setup_link()
  on.exit({
    for (k in names(old)) {
      if (is.na(old[[k]])) Sys.unsetenv(k)
      else do.call(Sys.setenv, setNames(list(old[[k]]), k))
    }
  }, add = TRUE)

  built <- tryCatch({
    Rcpp::sourceCpp(test_path("scm_gradient_driver.cpp"))
    Rcpp::sourceCpp(test_path("ad_certificate.cpp"))
    TRUE
  }, error = function(e) { message("build failed: ", conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  b0 <- 0.059; life <- 40

  ent <- k93_entry_gradient(c(1L, 6L, 7L), b_0 = b0, birth_rate = 20,
                            max_patch_lifetime = life, metric = 0L)

  # Reoptimising per-leaf certificate over the same schedule source (the certificate
  # driver refines from R). Its AD is the certified reference; its FD is a
  # reoptimising central difference on the resolved schedule.
  p0 <- scm_base_parameters("K93"); p0$max_patch_lifetime <- life
  pr <- add_strategies(p0, trait_matrix(b0, "b_0"), birth_rate = 20)
  base <- run_scm(pr, refine_schedule = TRUE)
  cert <- k93_allfield(base$parameters$node_schedule_times, base$parameters$ode_times,
                       b_0 = b0, birth_rate = 20, max_patch_lifetime = life, metric = 0L)

  nm  <- cert$names
  idx <- match(c("b_0", "d_0", "d_1"), nm)

  expect_equal(ent$value, ent$value_double, tolerance = 1e-10)
  # The entry's reverse gradient matches the certificate's certified reverse AD (same
  # schedule, same adjoint) ...
  expect_equal(ent$grad, cert$ad[idx], tolerance = 1e-6)
  # ... and the reoptimising FD (the correctness reference).
  expect_equal(ent$grad, cert$fd[idx], tolerance = 1e-3)
})
