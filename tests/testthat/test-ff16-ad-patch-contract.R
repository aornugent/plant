# AD-2 (#6): the FF16 Patch satisfies odelia's differentiable-System contract and
# closes the two seeding traps. Forward-mode (tape-free) check that seeding a
# trait through ad_parameters(), re-running prepare_strategy(), resetting and
# introducing a cohort carries an exact, nonzero derivative into a frozen derived
# quantity (height_0), and that the seed survives reset().
#
# Built out-of-tree via sourceCpp: the active FF16 physiology links against the
# freshly-built plant.so (which carries the forward-mode instantiation), and
# plant.so's XAD Tape symbols resolve against odelia.so.

find_plant_include <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  cand <- cand[file.exists(file.path(cand, "plant/ad_value.h"))]
  if (length(cand)) cand[1] else ""
}

find_plant_so <- function() {
  loaded <- getLoadedDLLs()
  if ("plant" %in% names(loaded)) {
    p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
    if (is.character(p) && length(p) == 1 && file.exists(p)) return(p)
  }
  ""
}

testthat::test_that("active FF16 Patch: seed -> prepare -> reset -> introduce carries the trait derivative", {
  plant_inc  <- find_plant_include()
  plant_so   <- find_plant_so()
  odelia_inc <- system.file("include", package = "odelia")
  odelia_so  <- system.file("libs", "odelia.so", package = "odelia")
  bh_inc     <- system.file("include", package = "BH")

  testthat::skip_if(!nzchar(plant_inc), "plant AD headers not found on include path.")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found to link the active Patch.")
  testthat::skip_if(!nzchar(odelia_so) || !file.exists(odelia_so),
                    "odelia shared library not found for tape linking.")

  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc)),
                         paste0("-I", shQuote(bh_inc))),
    PKG_LIBS = paste(shQuote(normalizePath(plant_so)),
                     shQuote(normalizePath(odelia_so))))

  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      #include <Rcpp.h>
      #include <XAD/XAD.hpp>
      #include <plant/control.h>
      #include <plant/strategy.h>
      #include <plant/parameters.h>
      #include <plant/models/ff16_strategy.h>
      #include <plant/models/ff16_environment.h>
      #include <plant/individual.h>
      #include <plant/internals.h>
      #include <plant/node.h>
      #include <plant/species.h>
      #include <plant/patch.h>

      // [[Rcpp::export]]
      Rcpp::NumericVector ad2_patch_contract() {
        using namespace plant;
        using ad     = xad::fwd<double>::active_type;
        using StratD = FF16_Strategy_<double>;
        using StratA = FF16_Strategy_<ad>;
        using PatchD = Patch<StratD, FF16_Environment>;

        Control ctrl;
        FF16_Environment env;

        // Resident (double) Patch, one FF16 species, lifted to the active scalar
        // via the AD-2 rebind_from contract.
        Parameters<StratD, FF16_Environment> pd;
        pd.strategies.push_back(StratD());
        pd.validate();
        PatchD patch_d(pd, env, ctrl);
        auto patch = patch_d.rebind_from<ad>();

        // Seed lma (column 0) with a unit forward tangent.
        std::vector<ad*> pars = patch.ad_parameters();
        const double n_params = (double)pars.size();
        xad::derivative(*pars[0]) = 1.0;

        patch.ad_prepare();                                        // trap 1
        const double seed_prep = xad::derivative(*patch.ad_parameters()[0]);
        patch.reset();                                             // trap 2
        const double seed_reset = xad::derivative(*patch.ad_parameters()[0]);
        patch.introduce_new_node(0);

        auto node = patch.at_species(0).node_begin();
        ad h0 = node->individual.state(HEIGHT_INDEX);
        const double d_h0 = xad::derivative(h0);

        // FD reference on the resident model (uniroot converged tightly).
        const double lma0  = xad::value(StratD().pars.lma);
        const double delta = lma0 * 1e-4;
        auto height0_at = [&](double lma) -> double {
          StratD s; s.pars.lma = lma;
          s.control.offspring_production_tol = 1e-13;
          s.control.offspring_production_iterations = 2000;
          s.prepare_strategy();
          return s.initial_height();
        };
        const double d_h0_fd =
          (height0_at(lma0 + delta) - height0_at(lma0 - delta)) / (2.0 * delta);

        return Rcpp::NumericVector::create(
            Rcpp::Named("n_params")   = n_params,
            Rcpp::Named("d_h0_ad")    = d_h0,
            Rcpp::Named("d_h0_fd")    = d_h0_fd,
            Rcpp::Named("seed_prep")  = seed_prep,
            Rcpp::Named("seed_reset") = seed_reset);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("undefined symbol", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD symbols unavailable in this session's plant.so.")
    stop(res)
  }

  r <- ad2_patch_contract()

  # All 32 FF16_Pars trait fields are exposed as differentiable handles.
  testthat::expect_equal(unname(r[["n_params"]]), 32)
  # The seeded trait's derivative reaches height_0 (nonzero) and matches FD.
  testthat::expect_gt(abs(r[["d_h0_ad"]]), 1e-6)
  testthat::expect_equal(unname(r[["d_h0_ad"]]), unname(r[["d_h0_fd"]]), tolerance = 1e-4)
  # The seed survives prepare and reset (traps 1 and 2).
  testthat::expect_equal(unname(r[["seed_prep"]]), 1)
  testthat::expect_equal(unname(r[["seed_reset"]]), 1)
})
