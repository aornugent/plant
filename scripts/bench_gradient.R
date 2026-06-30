# Same-machine A/B timing harness for the reverse-mode AD trait-gradient path
# (#472 scope B). Companion to scripts/bench_ab.R (SCM solve) and bench_tf24.R
# (TF24 hydraulics); this one times the GRADIENT calls.
#
# Why a bespoke harness: the gradient is two-pass -- pass 1 harvests a frozen
# schedule + per-RK-stage resident environment from a cached run_scm, pass 2 is
# a C++ AD replay + reverse sweep. Pass-1 harvest currently runs in R
# (ff16_harvest / tf24f_harvest, which rebuild the whole RcppR6 patch on each
# access -- O(stand size), see R/emergent_gradient.R:70-73). The refactor moves
# that harvest into C++, so this harness times the two passes SEPARATELY to make
# the win attributable:
#   run_ms     -- the resident run_scm(save_RK45_cache=TRUE) (context; timed once)
#   harvest_ms -- the R-side *_harvest() call alone           (the refactor target)
#   impl_ms    -- the C++ *_impl replay + sweep alone          (the AD kernel)
#   harvest_frac = harvest_ms / (harvest_ms + impl_ms)         (the headline)
#   val=<digest> -- a value digest of the Jacobian (same-session bit-identity)
#
# You cannot load two builds of `plant` in one R session, so build each and run
# this script once per build in a separate Rscript, interleaving to control for
# machine drift. Trust same-session ratios, not absolute numbers across sessions
# (see the profile-plant skill, .claude/skills/profile-plant/SKILL.md).
#
# Usage:
#   make compile                                   # build this branch
#   git worktree add -f /private/tmp/plant-other <ref> && ( cd /private/tmp/plant-other && make compile )
#   for r in 1 2; do
#     Rscript --no-init-file scripts/bench_gradient.R "$(pwd)"            "this-r$r"
#     Rscript --no-init-file scripts/bench_gradient.R /private/tmp/plant-other "other-r$r"
#   done
#
# Env vars:
#   BENCH_REPS   median over this many timed reps (default 7)
#   BENCH_CASES  comma-separated case subset (default all); e.g. "ff16_frozen,tf24f_census"
#   BENCH_SCALING set to 1 to also report impl scaling (1 vs all metrics; 4 vs all traits)

args  <- commandArgs(trailingOnly = TRUE)
path  <- if (length(args))     args[[1]] else "."
label <- if (length(args) > 1) args[[2]] else path
suppressMessages(pkgload::load_all(path, compile = FALSE, quiet = TRUE))

REPS    <- as.integer(Sys.getenv("BENCH_REPS", "7"))
SCALING <- nzchar(Sys.getenv("BENCH_SCALING", ""))

med_ms <- function(f, n = REPS) {
  ts <- numeric(n)
  for (i in seq_len(n)) ts[i] <- system.time(f())[["elapsed"]]
  median(ts) * 1000
}
# A cheap order-insensitive value digest: a Jacobian that changed anywhere moves it.
digest_jac <- function(x) {
  v <- if (is.list(x)) unlist(x[intersect(c("jacobian", "values", "d_state", "d_time",
                                            "time", "state", "gradient", "value"),
                                          names(x))], use.names = FALSE) else as.numeric(x)
  v <- v[is.finite(v)]
  sprintf("%.10g/%.10g/%d", sum(abs(v)), sum(v * seq_along(v)), length(v))
}

# ---- canonical cases (settings lifted verbatim from the test fixtures) --------

# FF16 frozen + offspring share one cached SCM (test-ff16-stand-gradient.R:7-12).
mk_ff16_basic <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# FF16 single-species resident coupled (test-ff16-stand-gradient.R:66-73).
mk_ff16_resident <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11)); p$max_patch_lifetime <- 60
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# FF16 two-species resident cross-species (test-ff16-stand-gradient.R:128-135).
mk_ff16_ms <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# TF24f census + resident share one cached SCM (test-tf24f-census-gradient.R:20-28).
mk_tf24f <- function(H = 4L, n = 9L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
}

# Each case: $run() builds the cached SCM (timed once -> run_ms); $harvest(ctx)
# the R-side harvest; $impl(ctx, h) the C++ kernel on a pre-harvested h; $public(ctx)
# the end-to-end public call (its return value is digested). harvest/impl NULL ->
# that sub-cost is reported NA (e.g. grow_individual's pass-1 is a bracket, not a harvest).
CASES <- list(
  ff16_frozen = list(
    run = mk_ff16_basic,
    mets = c("offspring_production", "LAI", "biomass", "size_moment"),
    tr = ff16_default_traits(),
    harvest = function(ctx) plant:::ff16_harvest(ctx, 1L, NULL),
    impl = function(ctx, h, mets, tr) plant:::ff16_stand_gradient_impl(
      h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr, mets, h$birth_rate,
      "frozen", list(), list(), h$patch_area, -1, -1),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr)),

  ff16_offspring = list(
    run = mk_ff16_basic,
    mets = "offspring_production",
    tr = ff16_default_traits(),
    harvest = function(ctx) plant:::ff16_harvest(ctx, 1L, NULL),
    impl = function(ctx, h, mets, tr) plant:::ff16_offspring_production_gradient_impl(
      h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr),
    public = function(ctx, mets, tr) offspring_production_gradient(ctx, traits = tr)),

  ff16_resident_coupled = list(
    run = mk_ff16_resident,
    mets = c("LAI", "biomass", "size_moment"),
    tr = c("a_p1", "lma", "a_l1"),
    harvest = function(ctx) plant:::ff16_harvest(ctx, 1L, NULL),
    impl = function(ctx, h, mets, tr) plant:::ff16_coupled_gradient_impl(
      h$pp, h$eh, h$sh, h$birth_step, h$ppsurv, h$ppsab, h$tw, tr, mets, h$birth_rate,
      h$nn_h, h$nn_c, h$patch_area),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr,
                                                    feedback = "resident")),

  ff16_resident_ms = list(
    run = mk_ff16_ms,
    mets = c("LAI", "size_moment"),
    tr = c("lma", "a_p1"),
    harvest = function(ctx) plant:::ff16_harvest_ms(ctx),
    impl = function(ctx, h, mets, tr) plant:::ff16_coupled_gradient_ms_impl(
      h$pp_list, h$eh, h$sh, h$birth_list, tr, mets, h$birth_rate, h$nn_h, h$nn_c,
      h$patch_area, 1L),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr,
                                                    species = 1L, feedback = "resident")),

  tf24f_census = list(
    run = mk_tf24f,
    mets = c("LAI", "biomass", "size_moment"),
    tr = c("vcmax_25", "lma", "a_l1", "K_s"),    # TF24f pars subset (test-tf24f-census-gradient.R:75)
    harvest = function(ctx) plant:::tf24f_harvest(ctx, 1L, NULL),
    impl = function(ctx, h, mets, tr) plant:::tf24f_census_gradient_ad_impl(
      h$pp, h$eh, h$sh, h$birth_step, h$birth_rate, h$k_acclim, h$use_ad_gradient,
      tr, mets, 1e-5),
    public = function(ctx, mets, tr) tf24f_census_gradient_ad(ctx, metrics = mets, traits = tr)),

  tf24f_resident = list(
    run = mk_tf24f,
    mets = c("LAI", "size_moment"),
    tr = c("vcmax_25", "lma", "a_l1", "K_s"),    # TF24f pars subset (test-tf24f-census-gradient.R:132)
    harvest = function(ctx) plant:::tf24f_harvest(ctx, 1L, NULL),
    impl = function(ctx, h, mets, tr) plant:::tf24f_coupled_gradient_impl(
      h$pp, h$eh, h$sh, h$birth_step, h$birth_rate, h$k_acclim, h$use_ad_gradient,
      tr, mets, h$nn_h, h$nn_c, h$patch_area, 1e-5),
    public = function(ctx, mets, tr) tf24f_resident_census_gradient_ad(ctx, metrics = mets,
                                                                       traits = tr)),

  grow_individual = list(
    run = function() {                     # the "run" is the individual + env, not an SCM
      list(indv = Individual("FF16", "FF16_Env")(FF16_Strategy()), env = Environment("FF16"),
           targets = c(2, 5, 10))
    },
    mets = NA, tr = ff16_default_traits(),
    harvest = NULL, impl = NULL,           # pass-1 is a bracket, not the resident harvest
    public = function(ctx, mets, tr) grow_individual_to_size_gradient(
      ctx$indv, ctx$targets, "height", ctx$env, traits = tr, time_max = 200))
)

# ---- run -----------------------------------------------------------------------

want <- Sys.getenv("BENCH_CASES", "")
sel  <- if (nzchar(want)) strsplit(want, ",")[[1]] else names(CASES)

for (nm in sel) {
  cs <- CASES[[nm]]
  if (is.null(cs)) { message("unknown case: ", nm); next }
  mets <- cs$mets; tr <- cs$tr

  run_ms <- med_ms(function() cs$run(), n = max(1L, REPS %/% 2L))   # solve is slow; fewer reps
  ctx <- cs$run()
  invisible(if (!is.null(cs$harvest)) cs$harvest(ctx))             # warm up
  invisible(cs$public(ctx, mets, tr))

  harvest_ms <- impl_ms <- NA_real_
  if (!is.null(cs$harvest)) {
    harvest_ms <- med_ms(function() cs$harvest(ctx))
    h <- cs$harvest(ctx)
    impl_ms <- med_ms(function() cs$impl(ctx, h, mets, tr))
  }
  public_ms <- med_ms(function() cs$public(ctx, mets, tr))
  val <- digest_jac(cs$public(ctx, mets, tr))
  hf  <- if (is.na(harvest_ms)) NA_real_ else harvest_ms / (harvest_ms + impl_ms)

  scaling <- ""
  if (SCALING && !is.null(cs$impl) && !identical(mets, NA) && length(mets) > 1L) {
    h <- cs$harvest(ctx)
    impl_1met  <- med_ms(function() cs$impl(ctx, h, mets[1], tr))
    impl_4tr   <- med_ms(function() cs$impl(ctx, h, mets, tr[seq_len(min(4L, length(tr)))]))
    scaling <- sprintf("|impl_1met_ms=%.3f|impl_4tr_ms=%.3f", impl_1met, impl_4tr)
  }

  cat(sprintf(
    "RESULT|%s|%s|run_ms=%.3f|harvest_ms=%.3f|impl_ms=%.3f|public_ms=%.3f|harvest_frac=%.3f|val=%s%s\n",
    label, nm, run_ms, harvest_ms, impl_ms, public_ms, hf, val, scaling))
}
