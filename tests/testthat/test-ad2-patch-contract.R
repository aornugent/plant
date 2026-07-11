# AD-2: the odelia System AD contract on FF16. The (double) model held in R is
# lifted to an active twin per gradient call (rebind_from), a chosen subset of
# one shared strategy's traits is seeded (ad_parameters, a fixed column order),
# and every cohort aliases that one strategy so the seed reaches late
# introductions. This checks the contract mechanics with forward-mode XAD
# (header-only, no reverse tape / odelia.so needed): value lift, the column
# order, pointer aliasing, forward-derivative seeding, and cohort sharing.
#
# The re-freeze of the *derived* quantities under the active scalar
# (prepare_strategy -> eta_c / area_leaf_0 carrying the trait derivative) lands
# with the FF16 invasion port, where the active assimilation path (QK
# active-bound integration + frozen-environment read) first compiles.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ff16_ad2_contract <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/models/ff16_strategy.h"))
  testthat::skip_if(!any(has_hdr),
                    "FF16 strategy header not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  testthat::skip_if(!nzchar(odelia_inc),
                    "odelia include (XAD) not found.")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc))))

  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <vector>
    #include <XAD/XAD.hpp>
    #include <plant/models/ff16_strategy.h>
    #include <plant/models/ff16_environment.h>
    #include <plant/individual.h>
    using fwd = xad::fwd<double>::active_type;

    // [[Rcpp::export]]
    int ff16_ad_n_parameters() {
      plant::FF16_Strategy_<fwd> sa;
      return static_cast<int>(sa.ad_parameters().size());
    }

    // [[Rcpp::export]]
    Rcpp::NumericVector ff16_ad_contract(double eta_in) {
      plant::FF16_Strategy sd;                 // double model held in R
      sd.pars.eta = eta_in;
      plant::FF16_Strategy_<fwd> sa = sd.rebind_from<fwd>();   // lift to active
      double lifted_eta = xad::value(sa.pars.eta);            // (1) lift round-trips

      std::vector<fwd*> ps = sa.ad_parameters();
      // eta is column 4 in FF16_Pars_::field_ptrs()
      xad::value(*ps[4]) = 42.0;                               // (2) alias into pars
      double aliased = xad::value(sa.pars.eta);
      xad::derivative(*ps[4]) = 1.0;                           // (3) seed derivative
      double seeded_deriv = xad::derivative(sa.pars.eta);

      // (4) a cohort built on the shared strategy sees the seed
      plant::FF16_Strategy_<fwd>::ptr strat =
          std::make_shared<plant::FF16_Strategy_<fwd>>(sa);
      plant::Individual<plant::FF16_Strategy_<fwd>, plant::FF16_Environment>
          cohort(strat);
      double cohort_eta = xad::value(strat->pars.eta);

      return Rcpp::NumericVector::create(
          Rcpp::_["lifted_eta"]   = lifted_eta,
          Rcpp::_["aliased_eta"]  = aliased,
          Rcpp::_["seeded_deriv"] = seeded_deriv,
          Rcpp::_["cohort_eta"]   = cohort_eta);
    }', verbose = FALSE)
}

testthat::test_that("FF16 exposes its trait fields in a fixed column order", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 AD-2 contract in pkgload load_all sessions.")
  compile_ff16_ad2_contract()

  # FF16_Pars_::field_ptrs() enumerates all 32 trait fields.
  expect_equal(ff16_ad_n_parameters(), 32L)
})

testthat::test_that("rebind_from lifts trait values and ad_parameters seeds them", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping FF16 AD-2 contract in pkgload load_all sessions.")
  compile_ff16_ad2_contract()

  r <- ff16_ad_contract(eta_in = 12.5)
  # rebind_from copies the double trait value onto the active scalar.
  expect_equal(unname(r["lifted_eta"]), 12.5)
  # ad_parameters() returns pointers into the one shared pars store.
  expect_equal(unname(r["aliased_eta"]), 42.0)
  # seeding a parameter's forward derivative is visible on that pars field.
  expect_equal(unname(r["seeded_deriv"]), 1.0)
  # a cohort built on the shared strategy sees the seeded value.
  expect_equal(unname(r["cohort_eta"]), 42.0)
})
