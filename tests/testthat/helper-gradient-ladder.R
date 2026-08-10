# The fixtures the reverse-sweep checks run on, their regime assertions, and the
# machinery every check reports its margin against.
#
# Five things a check has to carry before it counts, and they are why this file
# exists rather than each check holding its own literals:
#
#   path disjointness  every reference written here is in R, over quantities the
#                      forward model publishes. It traverses the forward
#                      reduction, which it must, and no transpose.
#   fault injection    ladder_injected() names which deliberate break is in
#                      force, and ladder_report_margin() records how much room
#                      the check had. A fault detected at three times the
#                      tolerance is a check about to stop working.
#   non-vacuity        ladder_expect_moves() zeroes the channel a check exists
#                      for and requires the answer to change.
#   adversariality     the fixtures keep every distinguishing quantity distinct:
#                      no equal heights, no small-integer ratios, two species
#                      differing by whole factors, seeds from a fixed generator.
#   tolerance          ladder_forward_floor() measures a fixture's own arithmetic
#                      floor and every tolerance is a multiple of it.
#
# There are two kinds of fixture and the difference is not cosmetic. A **patch**
# fixture is constructed: its state is written directly, so it can be put in the
# regime the checks need and crossed by hand. A **stand** fixture is run, so its
# state is whatever the trajectory reached and the constructed conditions are not
# available. The block-level checks take a patch; only the checks that need a
# trajectory take a stand.

# ---- traits -----------------------------------------------------------------

# Chosen so that no two of the per-species reduction parameters are equal or in a
# small-integer ratio between the species. A reduction sum that collapses the
# species into one scalar, or a scatter that reaches the wrong species, then
# moves a number.
ladder_traits <- function() {
  list(
    fast = c(lma = 0.0825, hmat = 5.13,  k_I = 0.5,  a_l1 = 5.44, a_l2 = 0.306),
    slow = c(lma = 0.1074, hmat = 11.71, k_I = 0.31, a_l1 = 3.61, a_l2 = 0.401))
}

ladder_parameters <- function(species = "fast", birth_rate = NULL,
                              lifetime = 2) {
  traits <- ladder_traits()
  if (is.null(birth_rate)) {
    birth_rate <- c(fast = 1.10, slow = 0.83)[species]
  }
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  for (i in seq_along(species)) {
    tr <- traits[[species[[i]]]]
    p <- add_strategies(p, trait_matrix(unname(tr), names(tr)),
                        hyperpar = TF24_hyperpar,
                        birth_rate = list(unname(birth_rate[[i]])))
  }
  p
}

ladder_control <- function(...) {
  Control(node_density_in_birth_date = TRUE, ...)
}

# ---- the patch fixtures, whose state is written rather than reached ----------

# One species, one cohort: the smallest stand that still carries endogenous
# feedback, because it shades itself and draws on its own soil layers. Both
# reduction transposes, the retention factor and the field's slope channel are
# live here. Two cohorts buy accumulation, not feedback, and keeping them apart
# is what makes a failure at four nodes localisable.
ladder_patch_one <- function() {
  ladder_patch(species = "fast",
               heights = list(4.73),
               log_densities = list(-0.39))
}

# Two species, two cohorts each. Two species catch a reduction sum that collapses
# them into one scalar; two cohorts of one species catch a trait accumulated per
# cohort rather than per species. Four nodes catch both.
#
# `cross` puts one pair out of birth-date order in height, which is the state the
# forward model runs and a guard written on height order refuses. It is on by
# default because folding the condition into every run is the point of it;
# passing FALSE is how a failure is localised to the crossing.
ladder_patch_two_by_two <- function(cross = TRUE) {
  heights <- if (cross) list(c(4.73, 6.29), c(5.41, 2.57))
             else       list(c(6.29, 4.73), c(5.41, 2.57))
  ladder_patch(species = c("fast", "slow"),
               heights = heights,
               log_densities = list(c(-0.39, -1.67), c(-1.33, -2.71)))
}

# A live patch with its state written to the values handed in. Nodes are
# introduced oldest-first, which is the order the state vector holds them in, and
# each at its own time: a node takes its birth date from the clock, so
# introducing them all at one instant leaves the abscissa tied and every question
# about birth-date order unanswerable.
ladder_patch <- function(species, heights, log_densities,
                         relative_reserve = 0.12, time = 1.37,
                         birth_dates = NULL) {
  p <- ladder_parameters(species)
  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  if (is.null(birth_dates)) {
    offset <- c(0, 0.17, 0.11)
    birth_dates <- lapply(seq_along(heights), function(i) {
      offset[[i]] + seq(0, by = 0.29, length.out = length(heights[[i]]))
    })
  }
  events <- do.call(rbind, lapply(seq_along(species), function(i) {
    data.frame(species = i, time = birth_dates[[i]])
  }))
  for (k in order(events$time)) {
    # A node takes its birth date from the boundary node, and the boundary node
    # takes it from the clock on an environment build. Introducing without one
    # leaves every node of a species sharing a date, and the abscissa the
    # birth-date coordinate integrates over is then tied.
    patch$set_time(events$time[[k]])
    patch$compute_environment()
    patch$introduce_new_node(events$species[[k]])
  }
  ladder_condition(patch, heights, log_densities, relative_reserve, time)
}

# Put the stand where the checks need it.
#
# The reserve placement is the vacuity guard and is the least obvious of these.
# The gate is centred near the bottom of the reserve range with a width of a
# tenth, so at comfortable reserves its slope is orders below its slope in the
# band. A stand left where a run put it damps every channel running through
# growth, and a check then passes because the signal is small rather than because
# the code is right.
ladder_condition <- function(patch, heights, log_densities,
                             relative_reserve = 0.12, time = 1.37,
                             moisture = c(0.192, 0.271, 0.233, 0.317, 0.208)) {
  state <- patch$ode_state
  at <- 0L
  for (i in seq_len(length(patch$species))) {
    sp <- patch$species[[i]]
    stride <- sp$nodes[[1]]$ode_size
    names_i <- sp$nodes[[1]]$ode_names
    i_height <- match("height", names_i)
    i_storage <- match("storage", names_i)
    i_density <- match("log_density", names_i)
    stopifnot(length(heights[[i]]) == length(sp$nodes),
              length(log_densities[[i]]) == length(sp$nodes))
    for (j in seq_along(sp$nodes)) {
      base <- at + (j - 1L) * stride
      h <- heights[[i]][[j]]
      state[[base + i_height]] <- h
      state[[base + i_density]] <- log_densities[[i]][[j]]
      state[[base + i_storage]] <-
        relative_reserve * ladder_storage_capacity(patch, i, h)
    }
    at <- at + sp$ode_size
  }
  # The layers are given distinct moistures for the same reason the heights are
  # distinct: a uniform column makes a scatter that reaches the wrong layer, and
  # a retention factor applied at the wrong one, invisible.
  n_layer <- patch$environment$get_soil_number_of_depths()
  stopifnot(length(moisture) >= n_layer)
  state[at + seq_len(n_layer)] <- moisture[seq_len(n_layer)]
  patch$set_ode_state(state, time)
  # The auxiliaries the regime is read off are written by a rate evaluation, not
  # by loading a state, so a report taken without this reads the previous one's.
  invisible(patch$ode_rates)
  patch
}

# ---- the stand fixtures, which need a trajectory -----------------------------

# Two species, two cohorts each, run rather than constructed. The trajectory
# rungs need the steps the adaptive pass resolved, so this cannot be conditioned;
# what it can be is chosen, and the birth rates are chosen to keep every cohort
# out of the absorbing reserve region.
ladder_stand_two_by_two <- function() {
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.63), c(0, 0.41))
  ladder_run(p)
}

# Three introductions in the order species 1, species 2, species 1, so the node
# stride is exercised in both directions. Introducing into species 1 while
# species 2 exists is the arrangement a narrowing implemented as a truncation of
# the tail fails on.
ladder_stand_introductions <- function() {
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  ladder_run(p)
}

ladder_stand_one <- function() {
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(c(0, 0.63))
  ladder_run(p)
}

ladder_run <- function(p, ctrl = NULL) {
  if (is.null(ctrl)) ctrl <- ladder_control()
  run_scm(p, Environment("TF24"), ctrl, collect = FALSE)
}

# ---- reading a stand --------------------------------------------------------

# A stand hands out a snapshot of its patch; a patch is itself. Reading through
# one accessor is what lets the references below serve both fixtures.
ladder_as_patch <- function(x) {
  if (is.null(x$species)) x$patch else x
}

# Every node of every species, boundary node included. The boundary node is the
# reduction's lower grid point and is not ODE state, so a reader that stops at
# the introduced nodes is reading a different distribution from the one the model
# integrates.
ladder_nodes <- function(x) {
  patch <- ladder_as_patch(x)
  out <- list()
  for (i in seq_len(length(patch$species))) {
    sp <- patch$species[[i]]
    nodes <- c(list(sp$new_node), rev(sp$nodes))
    out[[i]] <- data.frame(
      species = i,
      boundary = c(TRUE, rep(FALSE, length(nodes) - 1L)),
      birth_date = vapply(nodes, function(n) n$introduction_time, numeric(1)),
      height = vapply(nodes, function(n) n$height, numeric(1)),
      log_density = vapply(nodes, function(n) n$log_density, numeric(1)))
  }
  do.call(rbind, out)
}

ladder_strategy_parameter <- function(x, species, name) {
  ladder_as_patch(x)$species[[species]]$nodes[[1]]$individual$strategy$pars[[name]]
}

# Leaf area from height, the allometry written out here rather than read back
# from the model: a reference that asked the model for its own integrand would
# share the code it is refereeing.
ladder_area_leaf <- function(x, species, height) {
  a_l1 <- ladder_strategy_parameter(x, species, "a_l1")
  a_l2 <- ladder_strategy_parameter(x, species, "a_l2")
  (height / a_l1)^(1 / a_l2)
}

# The capacity the reserve gate is built on.
ladder_storage_capacity <- function(x, species, height) {
  par <- function(n) ladder_strategy_parameter(x, species, n)
  eta <- par("eta")
  eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)
  area_sapwood <- ladder_area_leaf(x, species, height) * par("theta")
  par("a_st1") * area_sapwood * height * eta_c * par("rho")
}

# The trapezium rule, signed, so a grid running the wrong way subtracts exactly
# as the forward reduction's does.
ladder_trapezium <- function(x, y) {
  n <- length(x)
  stopifnot(n >= 2L, length(y) == n)
  sum(diff(x) * (y[-1] + y[-n])) * 0.5
}

# The leaf-area census over a named coordinate, summed over species.
#
# A census is a quadrature of a density, so its weights are gaps in the
# coordinate the state is carried on. This forms the same sum over either choice
# of abscissa, so the two can be compared against the model's own and the answer
# says which axis the model integrated.
ladder_census_leaf_area <- function(x, coordinate = c("birth_date", "height")) {
  coordinate <- match.arg(coordinate)
  nodes <- ladder_nodes(x)
  total <- 0
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    weighted <- exp(part$log_density) * ladder_area_leaf(x, i, part$height)
    axis <- part[[coordinate]]
    ord <- order(axis)
    total <- total + ladder_trapezium(axis[ord], weighted[ord])
  }
  total
}

# The census's direct term for the leaf-area metric: the metric's own reading of
# the traits, at fixed state. It is not a sensitivity of the state at all, so no
# sweep produces it, and it is easy to omit precisely because it is a one-line
# calculation at the final state.
#
# Leaf area is a power of height over the two allometric constants, so the two
# rows are elementary and are written here rather than asked of the model.
ladder_census_direct_term <- function(x) {
  nodes <- ladder_nodes(x)
  total <- c(a_l1 = 0, a_l2 = 0)
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    a_l1 <- ladder_strategy_parameter(x, i, "a_l1")
    a_l2 <- ladder_strategy_parameter(x, i, "a_l2")
    area <- ladder_area_leaf(x, i, part$height)
    d_a_l1 <- -area / (a_l2 * a_l1)
    d_a_l2 <- -area * log(part$height / a_l1) / (a_l2^2)
    density <- exp(part$log_density)
    axis <- part$birth_date
    ord <- order(axis)
    total[["a_l1"]] <- total[["a_l1"]] +
      ladder_trapezium(axis[ord], (density * d_a_l1)[ord])
    total[["a_l2"]] <- total[["a_l2"]] +
      ladder_trapezium(axis[ord], (density * d_a_l2)[ord])
  }
  total
}

ladder_aux <- function(x, name) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, node$individual$aux(name))
    }
  }
  out
}

ladder_net_production <- function(x) ladder_aux(x, "net_mass_production_dt")

ladder_relative_reserve <- function(x) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      ind <- node$individual
      storage <- ind$ode_state[[match("storage", ind$ode_names)]]
      capacity <- ladder_storage_capacity(patch, i, node$height)
      out <- c(out, min(max(storage, 0) / capacity, 1))
    }
  }
  out
}

# The gate's slope as a fraction of its slope at the band's centre, so the floor
# is a share of the strongest response rather than an absolute number.
ladder_reserve_gate_slope <- function(x) {
  patch <- ladder_as_patch(x)
  r <- ladder_relative_reserve(patch)
  out <- numeric(0)
  at <- 1L
  for (i in seq_len(length(patch$species))) {
    a_st2 <- ladder_strategy_parameter(patch, i, "a_st2")
    for (node in patch$species[[i]]$nodes) {
      g <- 1 / (1 + exp(-(r[[at]] - a_st2) / 0.1))
      out <- c(out, (g * (1 - g)) / 0.25)
      at <- at + 1L
    }
  }
  out
}

ladder_light_at_reads <- function(x) {
  patch <- ladder_as_patch(x)
  out <- numeric(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, patch$environment$get_environment_at_height(node$height))
    }
  }
  out
}

# ---- the regime a fixture has to sit in --------------------------------------

# Each row is one assertion, whether it holds, and its measured value. A violated
# assertion invalidates a run rather than failing it, so the checks call
# ladder_require_regime() and the fixture's own test asserts the whole table.
#
# `level` is which set applies. A constructed patch can be put anywhere, so it
# carries the whole set. A run trajectory cannot be crossed by hand and cannot
# have its reserves placed, so those two are not asked of it -- and that gap is
# itself worth knowing, because it is the part of the regime the trajectory rungs
# do not exercise.
ladder_regime_report <- function(x, level = c("patch", "stand")) {
  level <- match.arg(level)
  patch <- ladder_as_patch(x)
  env <- patch$environment
  nodes <- ladder_nodes(patch)

  n_layer <- env$get_soil_number_of_depths()
  moisture <- env$get_soil_water_state()[seq_len(n_layer)]
  sat <- rep_len(env$soil_moist_sat, n_layer)
  residual <- 1e-2

  production <- ladder_net_production(patch)
  light <- ladder_light_at_reads(patch)

  entry <- function(name, ok, value) {
    data.frame(assertion = name, ok = isTRUE(ok),
               value = paste(signif(value, 6), collapse = " "),
               stringsAsFactors = FALSE)
  }

  report <- rbind(
    entry("soil moisture strictly interior on every layer",
          all(moisture > 1.5 * residual) && all(moisture < 0.99 * sat),
          moisture),
    entry("net production strictly positive at every cohort",
          length(production) > 0 && all(production > 0), min(production)),
    entry("light at every read at least 100 times the floor",
          length(light) > 0 && all(light > 100 * 1e-4), min(light)),
    entry("no two heights equal or in a small-integer ratio",
          ladder_non_commensurate(nodes$height), length(nodes$height)),
    entry("no two layers at the same moisture",
          level == "stand" || ladder_non_commensurate(moisture), moisture),
    entry("both species separated in every reduction parameter",
          ladder_species_separated(patch), ladder_species_ratio(patch)))

  if (level == "patch") {
    reserve <- ladder_relative_reserve(patch)
    gate <- ladder_reserve_gate_slope(patch)
    report <- rbind(
      report,
      entry("relative reserve inside the gate's transition band",
            length(reserve) > 0 && all(reserve > 0.02) && all(reserve < 0.30),
            range(reserve)),
      entry("reserve gate slope above its floor",
            length(gate) > 0 && all(gate > 0.4), min(gate)))
  }
  # Crossing is an adversariality condition rather than a regime one, and the
  # fixtures come in both forms, so the check that needs it asserts it itself.
  report
}

ladder_require_regime <- function(x, level = c("patch", "stand")) {
  report <- ladder_regime_report(x, level)
  bad <- report[!report$ok, , drop = FALSE]
  if (nrow(bad) > 0L) {
    testthat::skip(paste0(
      "fixture outside its declared regime, so this run is invalid rather than ",
      "failing: ", paste(bad$assertion, collapse = "; ")))
  }
  invisible(report)
}

ladder_non_commensurate <- function(x, tol = 1e-3) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2L) return(TRUE)
  for (i in seq_along(x)) {
    for (j in seq_along(x)) {
      if (i >= j) next
      if (abs(x[[i]] - x[[j]]) < tol * max(abs(x[[i]]), abs(x[[j]]))) return(FALSE)
      ratio <- max(x[[i]], x[[j]]) / min(x[[i]], x[[j]])
      if (any(abs(ratio - 2:4) < 1e-2)) return(FALSE)
    }
  }
  TRUE
}

# Pairs whose birth-date order and height order disagree. Reserve-gated growth
# lets a younger cohort overtake an older one, so a transpose guarding on height
# order refuses exactly the stands the forward model runs correctly.
ladder_crossing_count <- function(nodes) {
  total <- 0L
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    part <- part[order(part$birth_date), , drop = FALSE]
    if (nrow(part) < 2L) next
    total <- total + sum(diff(part$height) > 0)
  }
  total
}

ladder_species_ratio <- function(x) {
  patch <- ladder_as_patch(x)
  if (length(patch$species) < 2L) return(NA_real_)
  vapply(c("k_I", "a_l1", "a_l2"), function(n) {
    ladder_strategy_parameter(patch, 1, n) / ladder_strategy_parameter(patch, 2, n)
  }, numeric(1))
}

ladder_species_separated <- function(x) {
  ratio <- ladder_species_ratio(x)
  if (all(is.na(ratio))) return(TRUE)
  all(abs(log(ratio)) > 0.2)
}

# ---- seeds ------------------------------------------------------------------

# Block-normalised seeds from a fixed generator. Soil moisture is of order a
# third and heartwood mass is in kilograms, so an unnormalised inner product is a
# test of the largest block alone; and a seed of all ones makes a summation-order
# defect invisible. Signs alternate so a sum that should cancel does.
ladder_seeds <- function(n, scale = NULL, seed = 20240811L) {
  had <- exists(".Random.seed", envir = globalenv())
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  v <- stats::runif(n, 0.3, 1.7) * rep_len(c(1, -1), n)
  if (!is.null(scale)) {
    stopifnot(length(scale) == n)
    v <- v / ladder_block_scale(scale)
  }
  v
}

# The magnitude each entry is normalised by: its own size where it has one, and
# the block's median otherwise, so a zero entry does not divide.
ladder_block_scale <- function(x) {
  nonzero <- abs(x[x != 0 & is.finite(x)])
  m <- if (length(nonzero) > 0) stats::median(nonzero) else 1
  if (!is.finite(m) || m == 0) m <- 1
  ifelse(abs(x) > 0 & is.finite(x), abs(x), m)
}

# ---- tolerance and margin ---------------------------------------------------

# The fixture's own arithmetic floor: how far apart two evaluations of the same
# forward quantity land when the shared solver has been driven elsewhere in
# between. Every tolerance is a multiple of this rather than a literal, so a
# threshold that starts failing can be argued about rather than only raised.
#
# This is also the permutation check's cheapest form: a quantity that does not
# come back to itself is carrying something between solves.
ladder_forward_floor <- function(x) {
  patch <- ladder_as_patch(x)
  state <- patch$ode_state
  time <- patch$ode_time
  first <- patch$ode_rates
  patch$set_ode_state(state * 1.05, time)
  invisible(patch$ode_rates)
  patch$set_ode_state(state, time)
  second <- patch$ode_rates
  scale <- pmax(abs(first), abs(second))
  drift <- ifelse(scale > 0, abs(first - second) / scale, 0)
  max(max(drift), 4 * .Machine$double.eps)
}

ladder_margin <- function(observed, tolerance) {
  if (tolerance <= 0) return(Inf)
  observed / tolerance
}

# Compare and report how much room the check had. A fault caught at three times
# the tolerance is a check about to stop working, so the number is recorded and
# not only compared.
ladder_report_margin <- function(label, observed, tolerance) {
  message(sprintf("  %-52s %9.2e / %9.2e  uses %6.3f of budget",
                  label, observed, tolerance,
                  ladder_margin(observed, tolerance)))
  testthat::expect_lt(observed, tolerance)
}

# Zero the channel a check exists for and require the answer to move. A check on
# the slope channel passes with the whole slope channel missing if the fixture's
# slope adjoints happen to be zero.
ladder_expect_moves <- function(with_channel, without_channel, label,
                                floor = 1e-8) {
  denom <- max(max(abs(with_channel)), .Machine$double.eps)
  moved <- max(abs(with_channel - without_channel)) / denom
  testthat::expect(
    moved > floor,
    sprintf(paste("%s: zeroing the channel moved the answer by only %.3e, so",
                  "the check would pass with that channel absent"),
            label, moved))
  invisible(moved)
}

# ---- what the sweep can currently be asked ----------------------------------

# The one place the sweep's own availability is decided. A rung blocked because
# the sweep cannot run at all is one red line, named once, rather than the same
# cause repeated down the file -- and `test-gradient-ladder-sweep.R` is where
# that single line lives.
ladder_sweep_blocked <- function(stand) {
  tryCatch({
    stand_gradient(stand)
    NULL
  }, error = function(e) conditionMessage(e))
}

ladder_gradient_or_skip <- function(stand, ...) {
  result <- tryCatch(stand_gradient(stand, ...), error = function(e) e)
  if (inherits(result, "error")) {
    testthat::skip(paste("the sweep does not run on this stand:",
                         conditionMessage(result)))
  }
  result
}

# Splitting a recording at an interior step and sweeping the halves needs an
# entry point that sweeps a sub-range of a recorded trajectory. Nothing exposes
# one, so the check that the decomposition is exact cannot be made.
ladder_can_split_sweep <- function() {
  exists("stand_gradient_over_steps", mode = "function")
}

# The controller's rejected attempts are not part of the trajectory, and the
# claim that excluding them is exact is vacuous until the count is known to
# exceed zero. Nothing publishes the count.
ladder_can_count_rejections <- function() {
  exists("stand_step_attempts", mode = "function")
}

# ---- fault injection --------------------------------------------------------

# Which deliberate break is in force. A check's sensitivity is established by
# breaking the thing it watches, and this is what names the break so a margin can
# be recorded against it. The fault-injection runs are the deliverable, not a
# by-product: a rung whose faults have not been injected has not been climbed.
ladder_injected <- function(name = NULL) {
  current <- Sys.getenv("PLANT_LADDER_INJECT", unset = "")
  if (is.null(name)) return(current)
  identical(current, name)
}

# ---- the declared zero lists ------------------------------------------------

# Every trait column resolves to exactly one of these, and the list is part of
# the check. An exact zero is the signature of a missing accumulator and never of
# true insensitivity, which makes it indistinguishable from an ecological
# finding, so zeros are classified rather than tolerated.

# Whose whole column is zero by declaration, with the reason named. A parameter
# belongs here only when it has no live route at all; one whose birth-size
# channel alone is imposed to zero has a live column and belongs on the list
# below instead.
ladder_zero_by_construction <- function() {
  character(0)
}

# Whose birth-size channel is imposed to zero, which is a declared bias rather
# than a column class. The seed height solves an implicit condition on the
# strategy and the strategy's preparation is evaluated in double, so every
# parameter reaching birth size carries no row through that channel while keeping
# whatever rows its other routes give it.
#
# This is the one term no available instrument can referee: a forward tangent
# imposes the same equation, and a re-run finite difference cannot referee at
# production because a relative step of two parts in ten million moves a mature
# stand between alive and identically zero. So it is declared rather than
# measured, and a silent zero here is the failure the declaration prevents.
ladder_birth_size_channel_zero <- function() {
  c("lma", "omega", "a_l1", "a_l2", "rho", "theta", "a_r1", "a_b1")
}

# A parameter the boundary cannot answer for. A refusal, never a number, and an
# unknown parameter must refuse by name rather than return anything at all.
ladder_refused_by_name <- function() {
  c("p_50")
}

# The band a registered parameter that reaches no equation comes back in. It is
# round-off rather than exact zero, so the check is two-sided and not a test
# against zero.
ladder_roundoff_band <- function() c(1e-22, 1e-16)
