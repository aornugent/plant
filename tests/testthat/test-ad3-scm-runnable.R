# AD-3: the SCM is the runnable the odelia gradient driver accepts (Model A).
# Cohort introductions are interleaved between advance_fixed segments, so the
# SCM -- not the inner odelia Solver -- owns the replay; it satisfies the driver
# by adding get_system_ref(), a driver-managed tape, and rebind_from().
#
# This checks the runnable *surface* at compile time (header-only, no build of
# the active path, which lands at the FF16 invasion port). The behaviour of the
# double resident/mutant run through the fixed SCM::reset() is covered by
# test-scm.R / test-mutant.R.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ad3_scm_surface <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/scm.h"))
  testthat::skip_if(!any(has_hdr), "plant scm.h not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  testthat::skip_if(!nzchar(odelia_inc), "odelia include not found.")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc))))

  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <type_traits>
    #include <plant.h>
    using namespace plant;
    using scm_t = SCM<FF16_Strategy, FF16_Environment>;

    // get_system_ref() exposes the live patch the solver owns.
    static_assert(std::is_same_v<
        decltype(std::declval<scm_t&>().get_system_ref()),
        Patch<FF16_Strategy, FF16_Environment>&>);

    // rebind_from<S2>() is the double -> active mould (checked at S2=double, the
    // only scalar whose body compiles before the invasion port).
    static_assert(std::is_same_v<
        decltype(std::declval<const scm_t&>().template rebind_from<double>()),
        scm_t>);

    // The driver-managed tape keeps the SCM copyable (a shared_ptr), which the
    // RcppR6 bindings rely on.
    static_assert(std::is_copy_constructible_v<scm_t>);

    // [[Rcpp::export]]
    bool ad3_scm_runnable_surface_ok() { return true; }', verbose = FALSE)
}

testthat::test_that("SCM presents the odelia runnable surface", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping AD-3 SCM surface check in pkgload load_all sessions.")
  compile_ad3_scm_surface()
  expect_true(ad3_scm_runnable_surface_ok())
})
