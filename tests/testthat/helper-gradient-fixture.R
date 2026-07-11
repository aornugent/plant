# Shared machinery for the FF16 reverse-mode gradient regression baseline: the
# canonical resident run, its invasion + resident Jacobians, and the two-tier
# comparison against the recorded oracle. The AD test files and
# scripts/gradient_fixture.R build the Jacobians here so the recorded baseline and
# the tests that check it never drift apart.

# The compiled gradient entries live in the built plant.so and are namespace
# internal (RcppExports plus thin wrappers, not roxygen-exported). Under R CMD
# check the tests run inside the plant namespace so the bare names resolve; a
# pkgload::load_all session attaches an in-development DLL whose path carries
# "pkgload". The oracle is recorded against the installed build, so the AD test
# files defer to R CMD check and skip under load_all.
is_pkgload_dll_plant <- function() {
  loaded <- getLoadedDLLs()
  if (!("plant" %in% names(loaded))) return(FALSE)
  p <- tryCatch(loaded[["plant"]][["path"]], error = function(e) "")
  is.character(p) && length(p) == 1 && grepl("pkgload", p, fixed = TRUE)
}

# The single canonical case (matches the AD-6/AD-7/AD-8 gates): one FF16 resident,
# lma 0.1, 50-year patch, canopy recorded. Deterministic, so the Jacobians are a
# fixed function of the build.
gradient_fixture_traits <- c("lma", "hmat")
gradient_fixture_metrics <- list(
  invasion = c("offspring_production", "LAI", "biomass", "basal_area"),
  resident = c("LAI", "biomass", "basal_area"))

gradient_fixture_parameters <- function() {
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  add_strategies(p0, trait_matrix(0.1, "lma"), birth_rate = 1)
}

gradient_fixture_scm <- function() {
  ctrl <- Control()
  ctrl$save_RK45_cache <- TRUE
  run_scm(gradient_fixture_parameters(), Environment("FF16"), ctrl)
}

# The invasion and resident Jacobians the baseline records, computed the one way
# the tests and the generator both use.
gradient_fixture_jacobians <- function(scm = gradient_fixture_scm()) {
  list(
    invasion = invasion_gradient(scm, gradient_fixture_metrics$invasion,
                                 gradient_fixture_traits),
    resident = stand_gradient(scm, gradient_fixture_metrics$resident,
                              gradient_fixture_traits))
}

gradient_fixture_fingerprint <- function() {
  paste(R.version$platform, R.version$arch, sep = "/")
}

gradient_fixture_path <- function() {
  testthat::test_path("fixtures", "gradient-baseline.rds")
}

read_gradient_baseline <- function() {
  readRDS(gradient_fixture_path())
}

# Concatenate a Jacobian set into one numeric vector in a fixed order so fresh and
# recorded runs compare element-for-element (as.numeric on a matrix is column
# major, and the metric/trait order is fixed above).
flatten_gradient_jacobians <- function(j) {
  c(as.numeric(j$invasion$value), as.numeric(j$invasion$gradient),
    as.numeric(j$resident$value), as.numeric(j$resident$gradient))
}

# Cross-platform noise floor: floating-point reassociation on a different build
# moves the last few digits, so off the recording machine the oracle relaxes to
# the FD tolerance the AD gates already accept.
gradient_baseline_noise_floor <- 1e-4

# Two-tier match. On the machine that recorded the baseline the computation is
# deterministic, so a fresh run is bit-identical (the tight regression tier);
# elsewhere it must fall within the noise floor. Returns the verdict plus the
# observed spread so tests and the perturbation check can assert both directions.
gradient_values_agree <- function(fresh, recorded, fingerprint) {
  a <- as.numeric(fresh)
  b <- as.numeric(recorded)
  exact <- identical(fingerprint, gradient_fixture_fingerprint())
  if (length(a) != length(b)) {
    return(list(ok = FALSE, exact = exact, max_rel = Inf, n = length(b)))
  }
  max_rel <- max(abs(a - b) / pmax(abs(b), .Machine$double.eps))
  ok <- if (exact) all(a == b) else max_rel < gradient_baseline_noise_floor
  list(ok = ok, exact = exact, max_rel = max_rel, n = length(b))
}

expect_matches_gradient_baseline <- function(fresh, recorded, fingerprint, label) {
  r <- gradient_values_agree(fresh, recorded, fingerprint)
  testthat::expect_true(
    r$ok,
    info = sprintf("%s: %s tier, max rel diff %.3e over %d values", label,
                   if (r$exact) "bit-identity" else "noise-floor", r$max_rel, r$n))
}
