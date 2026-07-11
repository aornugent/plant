# AD-4: Patch satisfies odelia's Replayable concept, so the stepper drives the
# record hooks on the resident pass and the frozen-environment dispatch on the
# mutant/replay pass. The plant hooks were renamed to the concept's names
# (record_stage / record_ode_step / replay_step) and has_recorded_field() added;
# without the concept satisfied the frozen-read branch never fires and a mutant
# would silently recompute.
#
# This asserts the concept holds at compile time (also enforced in the compiled
# library via a static_assert in src/ff16_node.cpp). The behaviour of the double
# mutant replay through the renamed hooks is covered by test-mutant.R.

is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

compile_ad4_replayable <- function() {
  cand <- c(tryCatch(here::here("inst/include"), error = function(e) ""),
            system.file("include", package = "plant"))
  has_hdr <- file.exists(file.path(cand, "plant/patch.h"))
  testthat::skip_if(!any(has_hdr), "plant patch.h not found on include path.")
  plant_inc <- cand[has_hdr][1]
  odelia_inc <- system.file("include", package = "odelia")
  testthat::skip_if(!nzchar(odelia_inc), "odelia include not found.")
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(plant_inc)),
                         paste0("-I", shQuote(odelia_inc))))

  Rcpp::sourceCpp(code = '
    #include <Rcpp.h>
    #include <plant.h>

    static_assert(
        odelia::ode::Replayable<plant::Patch<plant::FF16_Strategy,
                                             plant::FF16_Environment>>,
        "Patch must satisfy odelia::ode::Replayable");

    // [[Rcpp::export]]
    bool ad4_patch_is_replayable() { return true; }', verbose = FALSE)
}

testthat::test_that("Patch satisfies odelia's Replayable concept", {
  testthat::skip_if(is_pkgload_dll_plant(),
    "Skipping AD-4 Replayable check in pkgload load_all sessions.")
  compile_ad4_replayable()
  expect_true(ad4_patch_is_replayable())
})
