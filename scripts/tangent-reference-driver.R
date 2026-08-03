# The whole-run tangent reference: seeds one trait, replays the double run's
# accepted ODE steps at the forward-mode scalar, and reports the census
# derivative against the adjoint and against a central difference of the same
# schedule. Writes the snapshot the three columns are read from later.
#
#   Rscript scripts/tangent-reference-driver.R <plant worktree> [out.csv] [traits]
#
# traits: comma-separated, default lma,hmat,k_I. A trait that fails is written
# into the snapshot with its error text rather than dropped.
#
# The library on .libPaths() must be the one the worktree was installed into,
# and odelia must be the library carrying odelia/hermite_interpolator.hpp,
# which plant/resource_spline.h includes. sourceCpp is given plant's installed
# .so because it cannot resolve plant's own symbols on its own.
#
# The lifetime is short on purpose. Forward mode measured over twenty times the
# double run, and at the production lifetime it did not finish in 2400 s of CPU
# twice; the configuration below runs in seconds and is where the reference is
# both valid and affordable.
args <- commandArgs(TRUE)
tree <- args[1]
if (is.na(tree)) stop("give the plant worktree as the first argument")
out_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else
  file.path(tree, "scripts", "tangent-reference.csv")
traits <- if (length(args) >= 3 && nzchar(args[3]))
  strsplit(args[3], ",")[[1]] else c("lma", "hmat", "k_I")

ODELIA_LIB <- "/home/user/lib-wave5/odelia"

MAX_PATCH_LIFETIME <- 2
INTRODUCTION_TIMES <- seq(0, 1, length.out = 8)
FD_RELATIVE_STEP <- 1e-5

Sys.setenv(CXX_STD = "CXX20")
Sys.setenv(PKG_CPPFLAGS = paste0(
  "-std=gnu++20 -DNDEBUG",
  " -isystem", system.file("include", package = "plant"),
  " -I", file.path(ODELIA_LIB, "include"),
  " -I", system.file("include", package = "Rcpp"),
  " -I", system.file("include", package = "BH")))
Sys.setenv(PKG_LIBS = paste("-Wl,--allow-shlib-undefined",
                            system.file("libs", "plant.so", package = "plant")))
suppressMessages({ library(plant); library(Rcpp) })
attach(asNamespace("plant"), name = "plant-internals")
sourceCpp(file.path(tree, "scripts", "tf24-tangent-census.cpp"))

say <- function(...) { cat(..., "\n", sep = ""); flush(stdout()) }

git_head <- function(dir) {
  out <- suppressWarnings(system2("git", c("-C", dir, "rev-parse", "HEAD"),
                                  stdout = TRUE, stderr = FALSE))
  if (length(out) == 1) out else NA_character_
}

# odelia is used as an installed library, not a checkout, so it has no commit to
# record. Its version and build stamp name the build and the digest over the
# headers is what a reader can re-derive to confirm they hold the same ones.
odelia_identity <- function(dir) {
  d <- read.dcf(file.path(dir, "DESCRIPTION"))
  headers <- sort(list.files(file.path(dir, "include"), recursive = TRUE,
                             full.names = TRUE))
  digest <- system2("sha256sum", headers, stdout = TRUE)
  sprintf("version %s; built %s; %d headers, digest %s",
          d[1, "Version"], d[1, "Built"], length(headers),
          substr(digest_of(paste(sub(" .*", "", digest), collapse = "")), 1, 16))
}

digest_of <- function(s) {
  f <- tempfile(); on.exit(unlink(f)); writeLines(s, f)
  sub(" .*", "", system2("sha256sum", f, stdout = TRUE))
}

# scm$run() and nothing else, so the introduction schedule is the one given
# above rather than an adaptively refined one.
build_scm <- function(trait_values = list()) {
  ctrl <- Control()
  s <- TF24_Strategy()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  if (length(trait_values)) {
    pars <- s$pars
    for (nm in names(trait_values)) pars[[nm]] <- trait_values[[nm]]
    s$pars <- pars
  }
  p <- scm_base_parameters("TF24")
  p$strategies <- list(s)
  p$max_patch_lifetime <- MAX_PATCH_LIFETIME
  p$node_schedule_times <- list(INTRODUCTION_TIMES)
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), ctrl)
  scm$run()
  list(scm = scm, ctrl = ctrl)
}

defaults <- TF24_Strategy()$pars

z <- build_scm()
scm <- z$scm

# Read the fields off the control the run was given rather than off the one the
# script meant to give it, so an edit above cannot leave the snapshot wrong.
non_default_control <- {
  d <- Control()
  used <- z$ctrl
  nm <- names(d)
  changed <- nm[!vapply(nm, function(k) identical(used[[k]], d[[k]]),
                        logical(1))]
  if (length(changed))
    paste(sprintf("%s=%s", changed, vapply(changed, function(k)
      paste(format(used[[k]]), collapse = ","), character(1))), collapse = " ")
  else "none"
}
nodes <- length(scm$patch$species[[1]]$nodes)
ode_size <- scm$patch$ode_size

# node_schedule_times[[1]] reads empty on a run whose schedule was set before
# it, so the configuration is confirmed from the resolved patch instead: the
# node count and ode_size are set by the introduction count, 73 = 1 + 8 * 9.
say(sprintf("plant %s", git_head(tree)))
say(sprintf("odelia %s\n       %s", ODELIA_LIB, odelia_identity(ODELIA_LIB)))
say(sprintf("lifetime %g  introductions %d  nodes %d  ode_size %d  time %.10f",
            MAX_PATCH_LIFETIME, length(INTRODUCTION_TIMES), nodes, ode_size,
            scm$patch$time))
say(sprintf("accepted ODE steps %d  Control non-default: %s",
            length(scm$ode_times), non_default_control))

base <- stand_census(scm)
metrics <- names(base)
say("census: ", paste(sprintf("%s %.15g", metrics, unlist(base)), collapse = "  "))

central_difference <- function(trait) {
  v0 <- defaults[[trait]]
  h <- FD_RELATIVE_STEP * v0
  up <- stand_census(build_scm(setNames(list(v0 + h), trait))$scm)
  down <- stand_census(build_scm(setNames(list(v0 - h), trait))$scm)
  (unlist(up) - unlist(down)) / (2 * h)
}

rows <- list()
for (trait in traits) {
  say("")
  say("--- ", trait, " ---")
  res <- try({
    adj <- stand_gradient(scm, traits = trait)$gradient[, trait]
    tt <- tangent_census_tf24(scm, trait, z$ctrl)
    fd <- central_difference(trait)
    list(adj = adj, tt = tt, fd = fd)
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    say("FAILED: ", conditionMessage(attr(res, "condition")))
    rows[[length(rows) + 1]] <- data.frame(
      trait = trait, metric = metrics, trait_value = defaults[[trait]],
      census_value = NA_real_, tangent = NA_real_, adjoint = NA_real_,
      central_fd = NA_real_, tangent_accepted_steps = NA_integer_,
      tangent_ode_size = NA_integer_,
      status = paste("failed:", conditionMessage(attr(res, "condition"))),
      stringsAsFactors = FALSE)
    next
  }

  tt <- res$tt
  # A tangent replay whose accepted step count or state size differs from the
  # double run's has not replayed that run, and its derivative belongs to a
  # different discretisation. Refuse it rather than record it.
  if (tt$accepted_steps != length(scm$ode_times) || tt$ode_size != ode_size)
    stop(sprintf("tangent replay diverged: steps %d vs %d, ode_size %d vs %d",
                 tt$accepted_steps, length(scm$ode_times), tt$ode_size, ode_size))
  # A trait the census does not respond to at this lifetime gates nothing: the
  # difference is identically zero and the two derivatives are at round-off, so
  # they agree for want of anything to disagree about. Say so on the row.
  dead <- all(res$fd == 0) &&
    max(abs(tt$tangent)) < 1e-12 * max(abs(unlist(base)))
  status <- if (dead)
    "no response: central difference identically zero, derivatives at round-off"
  else "ok"
  if (dead) say("no census response to this trait at lifetime ",
                MAX_PATCH_LIFETIME, "; this row gates nothing")

  value_rel <- max(abs(tt$value - unlist(base)) / pmax(abs(unlist(base)), 1e-300))
  say(sprintf("replay: steps %d  ode_size %d  introduced %d  time %.10f",
              tt$accepted_steps, tt$ode_size, tt$n_introduced, tt$time))
  say(sprintf("census value against the double run: rel %.3g", value_rel))
  for (i in seq_along(metrics)) {
    say(sprintf("%-18s tangent %22.15g  adjoint %22.15g  fd %20.12g  rel %10.4g",
                metrics[i], tt$tangent[i], res$adj[[metrics[i]]],
                res$fd[[metrics[i]]],
                abs(res$adj[[metrics[i]]] - tt$tangent[i]) /
                  max(abs(res$adj[[metrics[i]]]), abs(tt$tangent[i]))))
  }
  rows[[length(rows) + 1]] <- data.frame(
    trait = trait, metric = metrics, trait_value = defaults[[trait]],
    census_value = tt$value, tangent = tt$tangent,
    adjoint = unname(res$adj[metrics]), central_fd = unname(res$fd[metrics]),
    tangent_accepted_steps = tt$accepted_steps,
    tangent_ode_size = tt$ode_size, status = status, stringsAsFactors = FALSE)
}

snapshot <- do.call(rbind, rows)
header <- c(
  "# Whole-run tangent reference for the TF24 census, with the adjoint and a",
  "# central difference measured on the same run in the same session.",
  "#",
  "# The tangent and the adjoint share the leaf's supplied partials and both",
  "# drop d(height_0)/d(trait), so agreement between them is not evidence over",
  "# either. The central difference carries height_0 and neither derivative",
  "# does, so on the columns that otherwise agree it stands off by however much",
  "# that channel is worth for the trait: 2.6% for lma, 0.09% for k_I.",
  "#",
  sprintf("# date                  %s", format(Sys.Date())),
  sprintf("# plant commit          %s", git_head(tree)),
  sprintf("# plant worktree        %s", normalizePath(tree)),
  sprintf("# odelia library        %s", ODELIA_LIB),
  sprintf("# odelia identity       %s", odelia_identity(ODELIA_LIB)),
  sprintf("# R library             %s", .libPaths()[1]),
  sprintf("# max_patch_lifetime    %g", MAX_PATCH_LIFETIME),
  sprintf("# node_schedule_times   %s",
          paste(sprintf("%.15g", INTRODUCTION_TIMES), collapse = " ")),
  sprintf("# introductions         %d", length(INTRODUCTION_TIMES)),
  sprintf("# nodes                 %d", nodes),
  sprintf("# ode_size              %d", ode_size),
  sprintf("# accepted ODE steps    %d", length(scm$ode_times)),
  sprintf("# refine_schedule       FALSE (scm$run() only)"),
  sprintf("# Control               defaults, non-default fields: %s",
          non_default_control),
  sprintf("# birth_rate_y          1, is_variable_birth_rate FALSE"),
  sprintf("# traits                %s", paste(traits, collapse = " ")),
  sprintf("# central difference    relative step %g, re-run both sides",
          FD_RELATIVE_STEP))
writeLines(header, out_path)
suppressWarnings(write.table(snapshot, out_path, sep = ",", row.names = FALSE,
                             qmethod = "double", append = TRUE))
say("")
say("wrote ", out_path)
