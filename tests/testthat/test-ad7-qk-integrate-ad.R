# AD-7: qk::integrate_ad -- active-bound Gauss-Kronrod. The census height
# integral and the FF16 crown integral run from a fixed lower bound to the
# plant's own (active) height, so the upper bound carries a trait derivative;
# a frozen-node replay would miss the moving-bound sensitivity. This checks the
# rule is faithful to the double integrate() and that the bound derivative is
# correct (Leibniz: d/db of int_a^b f = f(b), plus central FD of the value).

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_qk_integrate_ad <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/qk.h"))
  testthat::skip_if(!any(has_hdr), "plant qk.h not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  plant_so <- system.file("libs", "plant.so", package = "plant")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found (QK rule tables).")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc))),
    PKG_LIBS = shQuote(normalizePath(plant_so)))

  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <XAD/XAD.hpp>
    #include <plant/qk.h>
    using plant::quadrature::QK;
    using fwd = xad::fwd<double>::active_type;

    // int_a^b x^2 dx via the double rule (reference).
    // [[Rcpp::export]]
    double qk_double_sq(double a, double b) {
      QK q(21);
      return q.integrate([](double x){ return x * x; }, a, b);
    }
    // Same integral via the active-bound rule (value).
    // [[Rcpp::export]]
    double qk_ad_sq_value(double a, double b) {
      QK q(21);
      return q.integrate_ad([](double x){ return x * x; }, a, b);
    }
    // d/db int_a^b x^2 dx via forward-mode on the upper bound.
    // [[Rcpp::export]]
    double qk_ad_sq_dvalue_db(double a, double b) {
      QK q(21);
      fwd bb(b); xad::derivative(bb) = 1.0;
      fwd r = q.integrate_ad([](fwd x){ return x * x; }, fwd(a), bb);
      return xad::derivative(r);
    }', verbose = FALSE)
}

testthat::test_that("integrate_ad matches double integrate() and is exact for polynomials", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping AD-7 integrate_ad in pkgload load_all sessions.")
  compile_qk_integrate_ad()

  a <- 0.3; b <- 2.7
  analytic <- (b^3 - a^3) / 3
  expect_equal(qk_ad_sq_value(a, b), qk_double_sq(a, b))   # bit-identical rule
  expect_equal(qk_ad_sq_value(a, b), analytic, tolerance = 1e-12)
})

testthat::test_that("integrate_ad carries the moving-bound derivative", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping AD-7 integrate_ad in pkgload load_all sessions.")
  compile_qk_integrate_ad()

  a <- 0.3; b <- 2.7
  # Leibniz: d/db int_a^b x^2 dx = f(b) = b^2.
  expect_equal(qk_ad_sq_dvalue_db(a, b), b^2, tolerance = 1e-10)
  # And consistent with a central finite difference of the value.
  h <- 1e-6
  fd <- (qk_ad_sq_value(a, b + h) - qk_ad_sq_value(a, b - h)) / (2 * h)
  expect_equal(qk_ad_sq_dvalue_db(a, b), fd, tolerance = 1e-6)
})
