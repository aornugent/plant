# AD-3 (#8): the active SCM is the runnable odelia's gradient driver calls.
# Forward-mode (tape-free) check that an active SCM<FF16_Strategy_<active>, ...>
# instantiates, exposes get_system_ref()/ad_parameters(), survives the reset fix
# (reset re-initialises from the seeded parameters instead of clobbering with the
# stored double snapshot), and replays a pinned (ode-times) schedule to completion
# at S=active with the seeded trait's derivative flowing through the stepped cohort
# state. Extends the AD-2 Patch check to a full SCM replay.
#
# Built out-of-tree via sourceCpp against the freshly-built plant.so (double
# symbols) + odelia.so (XAD tape linkage), exactly like the AD-2 contract test.

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

testthat::test_that("active SCM replays a pinned schedule to completion; seed survives reset", {
  plant_inc  <- find_plant_include()
  plant_so   <- find_plant_so()
  odelia_inc <- system.file("include", package = "odelia")
  odelia_so  <- system.file("libs", "odelia.so", package = "odelia")
  bh_inc     <- system.file("include", package = "BH")

  testthat::skip_if(!nzchar(plant_inc), "plant AD headers not found on include path.")
  testthat::skip_if(!nzchar(plant_so) || !file.exists(plant_so),
                    "plant shared library not found to link the active SCM.")
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
      #include <algorithm>
      #include <vector>
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
      #include <plant/scm.h>

      // [[Rcpp::export]]
      Rcpp::NumericVector ad3_scm_active_run() {
        using namespace plant;
        using ad     = xad::fwd<double>::active_type;
        using StratD = FF16_Strategy_<double>;
        using SCMD   = SCM<StratD, FF16_Environment>;

        Control ctrl;
        FF16_Environment env;

        Parameters<StratD, FF16_Environment> pd;
        pd.strategies.push_back(StratD());
        pd.max_patch_lifetime = 20.0;   // short-but-representative, several cohorts
        pd.validate();

        // Double resident run to obtain a schedule known to integrate stably (the
        // exact step times the adaptive stepper used); replay it fixed at S=active.
        SCMD scm_d(pd, env, ctrl);
        scm_d.run();
        std::vector<double> ode_times = scm_d.r_ode_times();
        std::sort(ode_times.begin(), ode_times.end());
        ode_times.erase(std::unique(ode_times.begin(), ode_times.end()),
                        ode_times.end());
        ode_times.front() = 0.0;
        ode_times.back()  = pd.max_patch_lifetime;

        // Lift to the active twin (rebind_from) and pin its schedule to the times.
        auto scm = scm_d.rebind_from<ad>();
        NodeSchedule ns = scm.r_node_schedule();
        ns.r_set_ode_times(ode_times);
        ns.r_set_use_ode_times(true);
        scm.r_set_node_schedule(ns);

        // Seed lma (column 0) with a unit forward tangent on the ACTIVE system,
        // reached through the runnable surface get_system_ref()/ad_parameters().
        std::vector<ad*> pars = scm.get_system_ref().ad_parameters();
        const double n_params = (double)pars.size();
        xad::derivative(*pars[0]) = 1.0;

        // The reset fix: re-init from the seeded parameters, do not clobber the
        // seed with the stored double snapshot.
        scm.reset();
        const double seed_reset =
          xad::derivative(*scm.get_system_ref().ad_parameters()[0]);

        // Full active replay to completion under the pinned schedule.
        scm.run();
        const double complete   = scm.complete() ? 1.0 : 0.0;
        const double final_time = scm.time();
        const double seed_end =
          xad::derivative(*scm.get_system_ref().ad_parameters()[0]);

        // The seeded trait derivative flowed through the stepped active cohort
        // state: the oldest cohort height carries a nonzero, finite derivative.
        const auto& sp = scm.get_system_ref().at_species(0);
        const double n_cohorts = (double)sp.size();
        double d_height = 0.0, height = 0.0;
        if (sp.size() > 0) {
          ad h = sp.node_begin()->individual.state(HEIGHT_INDEX);
          d_height = xad::derivative(h);
          height   = xad::value(h);
        }

        return Rcpp::NumericVector::create(
            Rcpp::Named("n_params")   = n_params,
            Rcpp::Named("complete")   = complete,
            Rcpp::Named("final_time") = final_time,
            Rcpp::Named("seed_reset") = seed_reset,
            Rcpp::Named("seed_end")   = seed_end,
            Rcpp::Named("n_cohorts")  = n_cohorts,
            Rcpp::Named("height")     = height,
            Rcpp::Named("d_height")   = d_height);
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    if (grepl("undefined symbol", conditionMessage(res), fixed = TRUE))
      testthat::skip("AD symbols unavailable in this session's plant.so.")
    stop(res)
  }

  r <- ad3_scm_active_run()

  # All 32 FF16_Pars trait fields are exposed through the runnable surface.
  testthat::expect_equal(unname(r[["n_params"]]), 32)
  # The active run replayed the pinned schedule to completion.
  testthat::expect_equal(unname(r[["complete"]]), 1)
  testthat::expect_equal(unname(r[["final_time"]]), 20, tolerance = 1e-8)
  testthat::expect_gt(unname(r[["n_cohorts"]]), 0)
  # The seed survives reset() and the whole replay (the reset fix).
  testthat::expect_equal(unname(r[["seed_reset"]]), 1)
  testthat::expect_equal(unname(r[["seed_end"]]), 1)
  # A grown cohort height is finite and carries the trait's derivative.
  testthat::expect_true(is.finite(r[["height"]]) && r[["height"]] > 0)
  testthat::expect_true(is.finite(r[["d_height"]]))
  testthat::expect_gt(abs(unname(r[["d_height"]])), 1e-8)
})
