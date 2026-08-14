## Light-field accuracy through EARLY stand development, and the stand it
## produces. Run once per build; every number is self-refereed against that
## build's own competition reduction, so no state crosses between versions.
source("/home/a/.claude/jobs/e02c60e6/tmp/d1-lib.R")

## The branch build carries a run-time placement switch; develop and the landed
## branch do not.
policy <- Sys.getenv("INTERP_POLICY")
if (nzchar(policy)) {
  args <- as.numeric(strsplit(Sys.getenv("INTERP_ARGS", "0.10,65,1,0"), ",")[[1]])
  interp_policy_set(policy, args[1], as.integer(args[2]), as.integer(args[3]), args[4])
}

label <- Sys.getenv("LABEL", "build")
eta <- FF16_Strategy()$pars$eta
ages <- as.numeric(strsplit(Sys.getenv("AGES", "0.5,1,2,3,5,8,12,20,40"), ",")[[1]])

cat(sprintf("### %s\n", label))
cat(sprintf("%6s %7s %8s %7s %12s %12s %14s %12s\n",
            "age", "canopy", "cohorts", "knots",
            "field max", "field med", "stem density", "offspring"))
for (a in ages) {
  scm <- tryCatch(ff16_scm(a, collect = FALSE), error = function(e) e)
  if (inherits(scm, "error")) {
    cat(sprintf("%6.1f   FAILED: %s\n", a, substr(conditionMessage(scm), 1, 50))); next
  }
  p <- scm$patch
  hs <- heights_of(p); ns <- densities_of(p)
  nk <- tryCatch({ v <- nrow(p$environment$light_availability$state)
                   if (is.null(v) || !length(v)) NA_integer_ else as.integer(v) },
                 error = function(e) NA_integer_)
  e <- field_error(p, eta)
  ## Stem density is the trapezium of n over the size distribution -- the
  ## quantity self-thinning is about.
  ok <- is.finite(hs) & is.finite(ns)
  o <- order(hs[ok])
  dens <- if (sum(ok) > 1) sum(diff(hs[ok][o]) * (head(ns[ok][o], -1) + tail(ns[ok][o], -1)) / 2)
          else NA_real_
  cat(sprintf("%6.1f %7.3f %8d %7s %12.3e %12.3e %14.6e %12.6f\n",
              a, max(hs), length(hs), if (is.na(nk)) "-" else as.character(nk),
              max(e$rel), median(e$rel), dens, scm$offspring_production))
}

## Is the field a function of the state, or of the build before it? Reach one
## state two ways and compare. develop rescales a carried knot set by default.
cat(sprintf("\n--- %s: is the field a function of the state alone? ---\n", label))
scm <- ff16_scm(10, collect = FALSE)
ran <- scm$patch
y <- ran$ode_state; t <- ran$ode_time
n <- vapply(ran$species, function(s) s$size, 0.0)
ran$set_ode_state(y, t)
after_run <- ran$environment$light_availability$state

fresh <- Patch("FF16", "FF16_Env")(scm$parameters, Environment("FF16"), Control())
ok <- tryCatch({
  fresh$set_state(t, y, n, numeric(0))
  fresh$set_ode_state(y, t)
  TRUE
}, error = function(e) { cat("  set_state unavailable:", substr(conditionMessage(e), 1, 60), "\n"); FALSE })
if (ok) {
  after_set <- fresh$environment$light_availability$state
  same_n <- identical(dim(after_run), dim(after_set))
  cat(sprintf("  knots run/set: %d / %d\n", nrow(after_run), nrow(after_set)))
  if (same_n) {
    cat(sprintf("  identical: %s   max |value diff|: %.3e\n",
                identical(after_run, after_set),
                max(abs(after_run[, 2] - after_set[, 2]))))
  } else {
    cat("  different knot counts: the field depends on how the state was reached\n")
  }
}
