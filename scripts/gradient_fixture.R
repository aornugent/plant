#!/usr/bin/env Rscript

# Record the validated FF16 invasion + resident gradient Jacobians to
# tests/testthat/fixtures/gradient-baseline.rds, the regression oracle every
# subsequent AD change is gated against. Run from the package root after a build:
#
#   Rscript scripts/gradient_fixture.R
#
# The Jacobians are deterministic, so the file is stable across re-runs on one
# build; the recorded platform fingerprint drives the two-tier tolerance in
# helper-gradient-fixture.R (bit-identity here, a noise floor elsewhere).

suppressMessages({
  if (requireNamespace("plant", quietly = TRUE) &&
      !isTRUE(getOption("plant.gradient_fixture.load_all"))) {
    library(plant)
  } else {
    pkgload::load_all(".", quiet = TRUE)
  }
})

source(file.path("tests", "testthat", "helper-gradient-fixture.R"))

jacobians <- gradient_fixture_jacobians()

baseline <- list(
  jacobians = jacobians,
  fingerprint = gradient_fixture_fingerprint(),
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  case = list(strategy = "FF16", lma = 0.1, birth_rate = 1,
              max_patch_lifetime = 50, save_RK45_cache = TRUE,
              traits = gradient_fixture_traits,
              metrics = gradient_fixture_metrics))

out <- file.path("tests", "testthat", "fixtures", "gradient-baseline.rds")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
saveRDS(baseline, out)

cat("wrote", out, "\n")
cat("fingerprint:", baseline$fingerprint, "\n")
cat("invasion gradient:\n"); print(jacobians$invasion$gradient)
cat("resident gradient:\n"); print(jacobians$resident$gradient)
