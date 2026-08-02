# V4, reference half: the re-run central finite difference of TF24's three
# census metrics and of R0, with respect to nine traits, on the identical
# resolved schedule. No adjoint is run here; this is the number the adjoint
# half of V4 is judged against.
#
#   Rscript scripts/v4-reference.R <plant worktree> [out.rds] [traits] [steps]
#
# traits: comma-separated subset of the nine, default all nine.
# steps:  comma-separated relative steps, default 1e-5 (the reference step).
#         Extra steps are the step-size check and are stored alongside.
#
# This is the reference half split out of scripts/v4-census-gradient.R, which
# keeps the adjoint comparison. The design (trait list, window, method) is that
# harness's; the differences are recorded in the header of the saved object.

args <- commandArgs(TRUE)
tree <- args[1]
if (is.na(tree)) stop("give the plant worktree as the first argument")
out_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else
  file.path(tree, "scripts", "v4-reference.rds")

suppressMessages(library(odelia))
pkgload::load_all(tree, quiet = TRUE)

MAX_PATCH_LIFETIME <- 105.32
LMA <- 0.1978791
ESTABLISHMENT_WINDOW <- c(3.222267, 8.544184)
FD_RELATIVE_STEP <- 1e-5

V4_TRAITS <- c("lma", "rho", "hmat", "theta", "a_l1", "k_I", "a_dG1", "K_s",
               "psi_crit")

traits <- if (length(args) >= 3 && nzchar(args[3]))
  strsplit(args[3], ",")[[1]] else V4_TRAITS
steps <- if (length(args) >= 4 && nzchar(args[4]))
  as.numeric(strsplit(args[4], ",")[[1]]) else FD_RELATIVE_STEP

say <- function(...) { cat(..., "\n", sep = ""); flush(stdout()) }

build_stand <- function() {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- MAX_PATCH_LIFETIME
  add_strategies(p, trait_matrix(LMA, "lma"))
}

# The same resolved schedule and the same ODE grid on every re-run, so the
# difference is of the trait and not of the discretisation.
run_on_schedule <- function(p, schedule_times, ode_times) {
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  scm$set_node_schedule_times(schedule_times)
  sched <- scm$node_schedule
  sched$ode_times <- ode_times
  sched$use_ode_times <- TRUE
  scm$node_schedule <- sched
  scm$run()
  scm
}

read_out <- function(scm) {
  c(stand_census(scm), R0 = scm$net_reproduction_ratios[[1]])
}

p <- build_stand()
t0 <- Sys.time()
base <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
base$run()
schedule <- base$node_schedule$all_times
ode_times <- base$ode_times
base_out <- read_out(base)
say(sprintf("base: %.1f s, %d accepted ODE steps, R0 = %.15g",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            length(ode_times), base_out[["R0"]]))

census_time <- base$time
in_window <- census_time >= ESTABLISHMENT_WINDOW[1] &&
  census_time <= ESTABLISHMENT_WINDOW[2]
say(sprintf("census at t = %.6f; establishment window [%.6f, %.6f]: %s",
            census_time, ESTABLISHMENT_WINDOW[1], ESTABLISHMENT_WINDOW[2],
            if (in_window) "INSIDE" else "outside"))
if (in_window) {
  stop("V4's census time is inside the establishment window; move it out ",
       "rather than mollifying the model")
}

sp <- base$patch$species[[1]]
nm <- sp$new_node$ode_names
state <- matrix(c(sp$ode_state, sp$new_node$ode_state), nrow = length(nm),
                dimnames = list(nm, NULL))
density <- exp(state["log_density", ])
n <- length(density)
closing_live <- density[[n]] > 0 && density[[n - 1]] > 0
say(sprintf("%d of %d cohorts at density exactly zero; boundary node %g: %s",
            sum(density == 0), n, density[[n]],
            if (closing_live) "closing trapezium is live"
            else "closing trapezium contributes nothing at this state"))

outputs <- names(base_out)
runs <- list()
fd <- array(NA_real_, c(length(outputs), length(traits), length(steps)),
            dimnames = list(outputs, traits, format(steps)))
one_sided <- array(NA_real_, c(length(outputs), length(traits), length(steps), 2),
                   dimnames = list(outputs, traits, format(steps),
                                   c("forward", "backward")))

for (tr in traits) {
  base_value <- p$strategies[[1]]$pars[[tr]]
  for (si in seq_along(steps)) {
    h <- steps[[si]] * max(abs(base_value), 1)
    vals <- list()
    for (sgn in c(-1, 1)) {
      pp <- p
      pp$strategies[[1]]$pars[[tr]] <- base_value + sgn * h
      tt <- Sys.time()
      scm <- run_on_schedule(pp, schedule, ode_times)
      got <- read_out(scm)
      ok_steps <- length(scm$ode_times) == length(ode_times)
      ok_sched <- isTRUE(all.equal(scm$node_schedule$all_times, schedule))
      runs[[length(runs) + 1]] <- list(
        trait = tr, step = steps[[si]], sign = sgn, value = got,
        n_ode_steps = length(scm$ode_times), steps_match = ok_steps,
        schedule_match = ok_sched,
        seconds = as.numeric(difftime(Sys.time(), tt, units = "secs")))
      say(sprintf("  %-9s h=%.0e %+d  %5.0f s  steps %d (%s) schedule %s",
                  tr, steps[[si]], sgn,
                  runs[[length(runs)]]$seconds, length(scm$ode_times),
                  if (ok_steps) "match" else "MISMATCH",
                  if (ok_sched) "match" else "MISMATCH"))
      vals[[as.character(sgn)]] <- got
      rm(scm); gc(FALSE)
    }
    fd[, tr, si] <- (vals[["1"]] - vals[["-1"]]) / (2 * h)
    one_sided[, tr, si, "forward"] <- (vals[["1"]] - base_out) / h
    one_sided[, tr, si, "backward"] <- (base_out - vals[["-1"]]) / h
    print(round(fd[, tr, si], 8))
    saveRDS(list(partial = TRUE, fd = fd, one_sided = one_sided, runs = runs),
            paste0(out_path, ".partial"))
  }
}

git_rev <- function(dir) {
  tryCatch(system2("git", c("-C", dir, "rev-parse", "HEAD"), stdout = TRUE),
           error = function(e) NA_character_)
}

ref <- list(
  what = "V4 reference: central finite difference of TF24 census metrics and R0",
  generated = Sys.time(),
  config = list(
    plant_commit = git_rev(tree),
    plant_worktree = normalizePath(tree),
    odelia_commit = "fdccd7b",
    odelia_lib = Sys.getenv("R_LIBS_USER"),
    build_flags = "CXX20FLAGS = -O2 -DNDEBUG -g0 (scripts/build/Makevars-O2)",
    model = list(strategy = "TF24", environment = "TF24_Env",
                 max_patch_lifetime = MAX_PATCH_LIFETIME,
                 lma = LMA, control = "Control() defaults",
                 refine_schedule = FALSE,
                 schedule = "base run's own schedule, pinned via
                   set_node_schedule_times + use_ode_times"),
    relative_steps = steps,
    reference_step = FD_RELATIVE_STEP,
    census_time = census_time,
    establishment_window = ESTABLISHMENT_WINDOW,
    traits = traits,
    outputs = outputs,
    excluded_traits = c("b", "c", "root_b", "root_c"),
    excluded_reason = "the cumulative vulnerability interpolant's knot count
      steps between 100 and 101 as these move by 1e-6 relative, so TF24's output
      is genuinely discontinuous in them and no finite difference is a reference"
  ),
  base_value = base_out,
  n_ode_steps = length(ode_times),
  cohorts = list(n = n, n_zero_density = sum(density == 0),
                 boundary_density = density[[n]],
                 closing_trapezium_live = closing_live),
  fd = fd,
  one_sided = one_sided,
  runs = runs)

saveRDS(ref, out_path)
say("wrote ", out_path)
