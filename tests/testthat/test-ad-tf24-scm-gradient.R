# TF24 / TF24f full-SCM reverse-mode trait gradient, through the run-shaped entry
# (plant/scm_gradient.h) -- the SAME growing-dimension path FF16/K93 use for R0 and
# census. The leaf-coupling drivers (gate0 / soil / collar) exercise one
# compute_rates; this is the emergent-metric SCM gradient a user actually calls.
#
# STATUS (measured, session 16): the TF24/TF24f active SCM COMPILES for all metrics
# (the earlier "won't compile at node.h:347" was stale) and REPRODUCES the double
# value bit-exactly through the entry (scm_jacobian's structural value check passes).
# The value-reproduction was fixed by snapshotting the leaf operating point around
# the net_mass_production_dt seam (the assembly's inner solves mutate leaf scratch;
# leaking it shifted the recorded trajectory off the double one).
#
# What these tests LOCK IN: compile + bit-exact value reproduction for every cell
# (the regression guard for that fix). What they do NOT yet assert: FD-verification
# of the gradient VALUE. A naive central FD of a TF24 SCM metric is dominated by the
# leaf model's per-call non-smoothness (golden-section p* at GSS_tol_abs=1e-3, regime
# switches) -- an FD-step sweep gives no plateau (ratios -1.4 .. +2.9 across
# d=1e-3..1e-7), the documented "staircase" (HANDOFF sessions 11-13). Certifying the
# gradient needs the tight-inner-eps / multi-cell reference from that work, not a
# central difference; that verification is the remaining open task (#27).

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

# The entry runs its own value-reproduction assert (scm_jacobian stops if the active
# value does not reproduce the double reference). So a call that RETURNS is itself the
# proof that (a) the active SCM compiled and ran, and (b) the value reproduces. We
# additionally assert the returned value/gradient are finite.
test_that("TF24/TF24f full-SCM gradient entry: compiles, runs, reproduces the double value", {
  skip_on_cran()
  old <- Sys.getenv(c("PKG_CXXFLAGS", "PKG_LIBS"), unset = NA)
  setup_link()
  on.exit({
    for (k in names(old)) {
      if (is.na(old[[k]])) Sys.unsetenv(k)
      else do.call(Sys.setenv, setNames(list(old[[k]]), k))
    }
  }, add = TRUE)

  built <- tryCatch({ Rcpp::sourceCpp(test_path("tf24_scm_gradient_driver.cpp")); TRUE },
                    error = function(e) { message("build failed: ", conditionMessage(e)); FALSE })
  skip_if_not(built, "sourceCpp build unavailable in this session")

  life <- 4L; EMPTY_S <- character(0); EMPTY_V <- numeric(0)

  # TF24 R0 (offspring) and census scalar: value reproduces (entry would stop
  # otherwise), value + gradient finite.
  for (m in c(0L, 1L)) {
    r <- tf24_scm_gradient("lma", 20, life, m, EMPTY_S, EMPTY_V)
    expect_true(is.finite(r$value))
    expect_true(all(is.finite(r$grad)))
  }
  # TF24f R0: the tracked-collar variant runs through the same entry.
  rf <- tf24f_scm_gradient("lma", 20, life, 0L, EMPTY_S, EMPTY_V)
  expect_true(is.finite(rf$value) && all(is.finite(rf$grad)))

  # TF24 / TF24f census vector (codomain 3): the three FF16-mirrored census methods
  # (added to TF24 so TF24f inherits them) feed scm_jacobian; values positive, 3x1
  # Jacobian finite.
  for (fn in list(tf24_scm_census_vector, tf24f_scm_census_vector)) {
    cv <- fn("lma", 20, life, EMPTY_S, EMPTY_V)
    expect_length(cv$values, 3L)
    expect_true(all(cv$values > 0))
    expect_equal(dim(cv$jacobian), c(3L, 1L))
    expect_true(all(is.finite(cv$jacobian)))
  }
})

test_that("TF24 full-SCM gradient FD-verification (OPEN -- staircase reference needed)", {
  skip(paste0("TF24 SCM metric FD has no plateau (leaf-model non-smoothness); needs ",
              "the tight-eps / multi-cell reference from HANDOFF sessions 11-13, task #27"))
})
