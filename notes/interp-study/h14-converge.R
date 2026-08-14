## What both families converge TO, so delta is chosen against a limit rather
## than against the incumbent. FF16's offspring production is the sharpest
## functional available: it moved 2e-03 between 65 and 257 canopy knots.
## Resolve siblings relative to this file, so the study runs from the repo.
STUDY <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
if (!file.exists(file.path(STUDY, "lib-field.R"))) STUDY <- "notes/interp-study"
source(file.path(STUDY, "lib-field.R"))

ff16_run <- function(lifetime = 20) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, lifetime * 0.9, length.out = 20))
  run_scm(p, use_ode_times = FALSE)
}

canopy_top <- function(res) {
  ## tallest cohort at the end of the run
  max(res$species[[1]]["height", , drop = TRUE], na.rm = TRUE)
}

cat("=== FF16 offspring production, both families refined ===\n")
cat(sprintf("%-22s %6s %20s %12s\n", "policy", "knots", "offspring", "vs finest"))
rows <- list()
for (s in list(list("canopy", 0.10,  65, 1, 0), list("canopy", 0.10, 129, 1, 0),
               list("canopy", 0.10, 257, 1, 0), list("canopy", 0.10, 513, 1, 0),
               list("canopy", 0.10,1025, 1, 0),
               list("fixed",  0.50,  65, 1, 0), list("fixed",  0.25,  65, 1, 0),
               list("fixed",  0.10,  65, 1, 0), list("fixed",  0.05,  65, 1, 0),
               list("fixed",  0.025, 65, 1, 0), list("fixed", 0.0125, 65, 1, 0))) {
  do.call(interp_policy_set, s)
  r <- try(ff16_run(), silent = TRUE)
  if (inherits(r, "try-error")) {
    cat(sprintf("%-22s   FAILED\n", sprintf("%s d=%.4f n=%d", s[[1]], s[[2]], s[[3]])))
    next
  }
  lab <- if (s[[1]] == "canopy") sprintf("canopy n=%d", s[[3]])
         else sprintf("fixed d=%.4f", s[[2]])
  nk <- length(r$patch$environment$light_availability$state[, 1])
  rows[[lab]] <- c(nk = nk, y = r$offspring_production)
  cat(sprintf("%-22s %6d %20.12e\n", lab, nk, r$offspring_production))
}

## the finest of each family is the reference for that family
fin_c <- rows[["canopy n=1025"]]
fin_f <- rows[["fixed d=0.0125"]]
cat(sprintf("\nfinest canopy (%d knots): %.12e\n", fin_c[["nk"]], fin_c[["y"]]))
cat(sprintf("finest fixed  (%d knots): %.12e\n", fin_f[["nk"]], fin_f[["y"]]))
cat(sprintf("the two families agree to: %.3e relative\n\n",
            abs(fin_f[["y"]] / fin_c[["y"]] - 1)))

cat(sprintf("%-22s %6s %12s\n", "policy", "knots", "vs limit"))
lim <- fin_f[["y"]]
for (lab in names(rows))
  cat(sprintf("%-22s %6d %12.3e\n", lab, rows[[lab]][["nk"]],
              rows[[lab]][["y"]] / lim - 1))
