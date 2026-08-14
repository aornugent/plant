## Every placement in ONE binary, so the only thing that differs is the grid.
## Reports light-field error through early stand development, and the thinning
## trajectory's deviation from a converged reference.
Sys.setenv(PLANT_BUILD = "/home/a/dev/plant-dev/plant/.claude/worktrees/interp-design")
source("/home/a/.claude/jobs/e02c60e6/tmp/d1-lib.R")

eta <- FF16_Strategy()$pars$eta
ages <- c(0.5, 1, 2, 3, 4, 5, 6, 8, 12, 20, 40)

settings <- list(
  `canopy x h_max, 65`  = list("canopy", 0.10,   65, 1, 0),
  `canopy x h_max, 257` = list("canopy", 0.10,  257, 1, 0),
  `fixed d=0.100`       = list("fixed",  0.100,  65, 1, 0),
  `fixed d=0.050`       = list("fixed",  0.050,  65, 1, 0),
  `fixed d=0.025`       = list("fixed",  0.025,  65, 1, 0),
  `fixed d=0.0125 (ref)`= list("fixed",  0.0125, 65, 1, 0))

## stand state at an age, under one placement
at_age <- function(cfg, a) {
  do.call(interp_policy_set, cfg)
  scm <- ff16_scm(a, collect = FALSE)
  p <- scm$patch
  hs <- heights_of(p); ns <- densities_of(p)
  ok <- is.finite(hs) & is.finite(ns) & hs > 0
  o <- order(hs[ok])
  dens <- sum(diff(hs[ok][o]) * (head(ns[ok][o], -1) + tail(ns[ok][o], -1)) / 2)
  e <- field_error(p, eta)
  list(knots = nrow(p$environment$light_availability$state),
       canopy = max(hs), n = sum(ok),
       field_max = max(e$rel), field_med = median(e$rel),
       density = dens, offspring = scm$offspring_production,
       A0 = p$compute_competition(0))
}

res <- list()
for (nm in names(settings))
  res[[nm]] <- lapply(ages, function(a) at_age(settings[[nm]], a))

cat("=== light-field error through stand development (max relative, crown mean) ===\n")
cat(sprintf("%6s %7s", "age", "canopy"))
for (nm in names(settings)) cat(sprintf(" %13s", substr(nm, 1, 13)))
cat("\n")
for (i in seq_along(ages)) {
  cat(sprintf("%6.1f %7.3f", ages[[i]], res[[1]][[i]]$canopy))
  for (nm in names(settings)) cat(sprintf(" %13.3e", res[[nm]][[i]]$field_max))
  cat("\n")
}

cat("\n=== knots carried ===\n")
cat(sprintf("%6s", "age"))
for (nm in names(settings)) cat(sprintf(" %13s", substr(nm, 1, 13)))
cat("\n")
for (i in seq_along(ages)) {
  cat(sprintf("%6.1f", ages[[i]]))
  for (nm in names(settings)) cat(sprintf(" %13d", res[[nm]][[i]]$knots))
  cat("\n")
}

ref <- "fixed d=0.0125 (ref)"
for (q in c("density", "A0", "offspring")) {
  cat(sprintf("\n=== %s: relative deviation from the converged reference ===\n", q))
  cat(sprintf("%6s %14s", "age", "reference"))
  for (nm in setdiff(names(settings), ref)) cat(sprintf(" %13s", substr(nm, 1, 13)))
  cat("\n")
  for (i in seq_along(ages)) {
    r <- res[[ref]][[i]][[q]]
    cat(sprintf("%6.1f %14.7e", ages[[i]], r))
    for (nm in setdiff(names(settings), ref)) {
      v <- res[[nm]][[i]][[q]]
      cat(sprintf(" %13.3e", if (abs(r) > 0) abs(v / r - 1) else NA_real_))
    }
    cat("\n")
  }
}
