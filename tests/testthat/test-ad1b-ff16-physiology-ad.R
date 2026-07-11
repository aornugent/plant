# AD-1b: the FF16 physiology compiles and differentiates at an active scalar --
# net_mass_production_dt through the templated allocation/respiration/turnover
# and the crown assimilation integral (qk::integrate_ad over the active height,
# reading the frozen light via the active-query spline eval). This is the
# strategy half of the active port; the demographic layer (Node/Species) and the
# end-to-end invasion FD gate land with AD-6.
#
# Reverse-mode, so it links odelia.so (the XAD tape) and plant.so (the double
# QK / Control symbols); the active FF16_Strategy_<adj> instantiation is compiled
# from the source .cpp, so the test runs against the source tree (dev/load_all),
# not the installed package.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_physiology_ad <- function() {
  src <- tryCatch(here::here("src"), error = function(e) "")
  inc <- tryCatch(here::here("inst/include"), error = function(e) "")
  testthat::skip_if(!nzchar(src) ||
                    !file.exists(file.path(src, "ff16_strategy.cpp")),
                    "FF16 strategy source not available (installed package).")
  odelia_inc <- system.file("include", package = "odelia")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  plant_so <- system.file("libs", "plant.so", package = "plant")
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found (QK / Control symbols).")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(inc)), paste0("-I", shQuote(odelia_inc))),
    PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                     shQuote(normalizePath(odelia_so))))

  Rcpp::sourceCpp(code = sprintf('
    #include <Rcpp.h>
    #include <XAD/XAD.hpp>
    #include <plant/models/ff16_strategy.h>
    #include <plant/models/ff16_environment.h>
    #include "%s/ff16_strategy.cpp"
    using adt = xad::adj<double>::active_type;

    // Reverse-mode d(net_mass_production_dt)/d(lma) at a fixed height and
    // (constant) light, through the templated FF16 physiology.
    // [[Rcpp::export]]
    double ff16_net_grad_lma(double lma0, double height, double light) {
      using ad = xad::adj<double>;
      ad::tape_type tape;
      plant::FF16_Strategy_<adt> s;
      adt lma = lma0; tape.registerInput(lma); tape.newRecording();
      s.pars.lma = lma; s.prepare_strategy();
      plant::FF16_Environment env; env.set_fixed_environment(light);
      adt h = height, al = s.area_leaf(h), hinv = adt(1.0) / h;
      adt net = s.net_mass_production_dt(env, h, al, hinv);
      tape.registerOutput(net); xad::derivative(net) = 1.0; tape.computeAdjoints();
      return xad::derivative(lma);
    }
    // Double-path value for the finite-difference reference.
    // [[Rcpp::export]]
    double ff16_net_value_lma(double lma0, double height, double light) {
      plant::FF16_Strategy_<double> s; s.pars.lma = lma0; s.prepare_strategy();
      plant::FF16_Environment env; env.set_fixed_environment(light);
      double h = height, al = s.area_leaf(h), hinv = 1.0 / h;
      return s.net_mass_production_dt(env, h, al, hinv);
    }', here::here("src")), verbose = FALSE)
}

testthat::test_that("active FF16 net_mass_production_dt gradient matches finite differences", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping AD-1b physiology AD in pkgload load_all sessions.")
  compile_ff16_physiology_ad()

  lma0 <- 0.1978791; height <- 5.0; light <- 0.8
  g <- ff16_net_grad_lma(lma0, height, light)
  h <- 1e-6 * lma0
  fd <- (ff16_net_value_lma(lma0 + h, height, light) -
         ff16_net_value_lma(lma0 - h, height, light)) / (2 * h)
  expect_equal(g, fd, tolerance = 1e-5)
})
