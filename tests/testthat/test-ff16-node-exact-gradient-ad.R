# Milestone C increment 10 (#472 scope B / #537 A1): wire the exact AD
# growth-rate gradient into the live solver path. Node::growth_rate_gradient,
# when control.node_gradient_exact_ad is set and the strategy provides an AD
# gradient, returns Individual::growth_rate_gradient_exact (the strategy's exact
# d(growth rate)/d(height)) instead of the finite difference; otherwise it falls
# back to FD (the default; non-FF16 strategies always fall back via the
# Strategy<E> base returning NA). Here we check, at the Individual level, that the
# exact gradient matches a fine FD of the REAL growth_rate_given_height (the
# actual compute_rates path), in a varying light profile.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_node_exact <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/models/ff16_strategy.h"))
  testthat::skip_if(!any(has_hdr), "FF16 headers not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  plant_so <- system.file("libs", "plant.so", package = "plant")
  odelia_so <- system.file("libs", "odelia.so", package = "odelia")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found.")
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")
  testthat::skip_if(!nzchar(system.file("include", package = "BH")),
                    "BH include dir not resolvable (e.g. R CMD check sandbox).")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc)),
                         paste0("-I", shQuote(system.file("include", package = "BH")))),
    PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                     shQuote(normalizePath(odelia_so))))

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      #include <vector>
      #include <plant/models/ff16_strategy.h>
      #include <plant/individual.h>
      typedef plant::Individual<plant::FF16_Strategy, plant::FF16_Environment> Ind;

      static plant::FF16_Environment varying_env() {
        plant::FF16_Environment env;
        std::vector<double> st = {0,5,10,15,22, 0.30,0.50,0.68,0.85,1.0};
        env.r_init_interpolators(st);
        return env;
      }

      // Exact AD gradient (the value Node returns when the flag is on) and a fine
      // FD of the REAL growth_rate_given_height (the compute_rates path Node FDs).
      // [[Rcpp::export]]
      Rcpp::NumericVector node_exact_vs_fd(double height) {
        plant::FF16_Strategy s; s.prepare_strategy();  // default deep-crown
        Ind ind = plant::make_individual<plant::FF16_Strategy,
                                         plant::FF16_Environment>(s);
        plant::FF16_Environment env = varying_env();
        ind.set_state(HEIGHT_INDEX, height);
        const double exact = ind.growth_rate_gradient_exact(env);
        const double e = 1e-6;
        const double fd = (ind.growth_rate_given_height(height + e, env) -
                           ind.growth_rate_given_height(height - e, env)) / (2 * e);
        return Rcpp::NumericVector::create(exact, fd);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("active_tape_", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD symbols unavailable in this load_all session.")
    stop(res)
  }
}

testthat::test_that("Node-path exact AD gradient matches FD of the real compute_rates", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 node exact-gradient AD in pkgload load_all sessions.")
  compile_ff16_node_exact()

  for (height in c(4, 8, 12)) {
    r <- node_exact_vs_fd(height)
    expect_equal(r[1], r[2], tolerance = 1e-6)
  }
})
