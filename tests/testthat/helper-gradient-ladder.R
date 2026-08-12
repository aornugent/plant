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
ladder_patch_two_by_two <- function(cross = TRUE, parameters = NULL) {
  heights <- if (cross) list(c(4.73, 6.29), c(5.41, 2.57))
             else       list(c(6.29, 4.73), c(5.41, 2.57))
  ladder_patch(species = c("fast", "slow"),
               heights = heights,
               log_densities = list(c(-0.39, -1.67), c(-1.33, -2.71)),
               parameters = parameters)
}

# A live patch with its state written to the values handed in. Nodes are
# introduced oldest-first, which is the order the state vector holds them in, and
# each at its own time: a node takes its birth date from the clock, so
# introducing them all at one instant leaves the abscissa tied and every question
# about birth-date order unanswerable.
ladder_patch <- function(species, heights, log_densities,
                         relative_reserve = 0.12, time = 1.37,
                         birth_dates = NULL, parameters = NULL) {
  p <- if (is.null(parameters)) ladder_parameters(species) else parameters
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

# The introductions fixture with every species' establishment constant scaled, so
# a recruit can be put where the establishment probability is stiff.
#
# The probability is P^2/(P^2 + k^2) times a decay in patch age, with
# k = a_d0 * a_0, so its derivative in P peaks at P = k/sqrt(3) -- ratio 1 below.
# At the shipped constant the fixture's recruits sit at ratio 36 and 21, which is
# the flat region where the check would prove nothing; the scale is what moves
# them onto the peak.
ladder_stand_marginal_recruit <- function(scale = 30) {
  p <- ladder_parameters(c("fast", "slow"))
  strategies <- p$strategies
  for (i in seq_along(strategies)) {
    pars <- strategies[[i]]$pars
    pars$a_d0 <- pars$a_d0 * scale
    strategies[[i]]$pars <- pars
  }
  p$strategies <- strategies
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  ladder_run(p)
}

# Where each species' recruit sits on the establishment probability's stiff band,
# as a multiple of the potential at which its derivative peaks. One is the peak.
ladder_recruit_stiffness <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)), function(i) {
    node <- patch$species[[i]]$new_node
    production <- node$individual$aux("net_mass_production_dt")
    k <- ladder_strategy_parameter(patch, i, "a_d0") *
      ladder_area_leaf(patch, i, node$height)
    production / (k / sqrt(3))
  }, numeric(1))
}

ladder_recruit_establishment <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)), function(i) {
    patch$species[[i]]$new_node$individual$establishment_probability(
      patch$environment)
  }, numeric(1))
}

# The smallest stand that shows the allometric-constant disagreement, and its
# control: one species, four tenths of a year, with and without a second cohort.
#
# It exists so work on that defect has a loop measured in seconds. The pair runs in
# about six seconds each against forty-three for the two-species fixture the
# measurements were first taken on, and it separates the defect by a factor of six
# thousand -- 1.4e-03 against 2.2e-07 on the trajectory term. A second cohort
# cannot be had without an introduction after t = 0, so the control differs by that
# introduction as well as by the cohort; every other check here has excluded the
# introduction, which is what makes the pair a usable control rather than a
# confound.
ladder_stand_allometric_probe <- function(two_cohorts = TRUE) {
  p <- ladder_parameters("fast", lifetime = 0.4)
  p$node_schedule_times <- list(if (two_cohorts) c(0, 0.29) else 0)
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

# Each species' inflow boundary density, as the patch currently holds it. The
# boundary node is not ODE state, so which evaluation of the condition a patch is
# carrying is a property of what was last done to it, and that is what the two
# readers below are for.
ladder_boundary_density <- function(x) {
  patch <- ladder_as_patch(x)
  vapply(seq_len(length(patch$species)),
         function(i) patch$species[[i]]$new_node$log_density, numeric(1))
}

# The inflow condition at each of its two evaluations, from one recorded state.
#
# A stage evaluates it twice. `in_field` is the one a state load produces, taken
# in the field with every species' boundary interval left off -- the field is then
# rebuilt including it, so this is the value the reductions were built on.
# `in_uptake` is the one a rate evaluation produces, taken in that rebuilt field;
# it is the value the water aggregation reads, the value an introduced node
# inherits, and the value a census reports.
#
# They are the same function at different arguments, so nothing about either number
# says which one a caller is holding.
ladder_boundary_evaluations <- function(patch, state, time) {
  patch$set_ode_state(state, time)
  in_field <- ladder_boundary_density(patch)
  invisible(patch$ode_rates)
  list(in_field = in_field, in_uptake = ladder_boundary_density(patch))
}

# A gradient's columns with their species index stripped. A check written about a
# strategy's parameter finds it in either species' column through this; a check
# about one species' column names it in full.
ladder_bare_traits <- function(x) sub("^[0-9]+\\.", "", x)

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
#
# One entry per species per parameter, named as the gradient's columns are. A
# species' own row reads only its own cohorts, so summing the species and
# comparing the total against a column would compare a sum of two things against
# one of them.
ladder_census_direct_term <- function(x) {
  nodes <- ladder_nodes(x)
  out <- numeric(0)
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
    term <- c(a_l1 = ladder_trapezium(axis[ord], (density * d_a_l1)[ord]),
              a_l2 = ladder_trapezium(axis[ord], (density * d_a_l2)[ord]))
    out <- c(out, stats::setNames(term, paste0(i, ".", names(term))))
  }
  out
}

# The rows of one rate evaluation, named. Nodes come species-major in the order
# the state vector holds them and the environment closes it, so a Jacobian row
# index can be read as the rate it belongs to.
ladder_rate_names <- function(x) {
  patch <- ladder_as_patch(x)
  out <- character(0)
  for (i in seq_len(length(patch$species))) {
    for (node in patch$species[[i]]$nodes) {
      out <- c(out, node$ode_names)
    }
  }
  n <- length(patch$ode_state)
  stopifnot(length(out) <= n)
  c(out, paste0("environment_", seq_len(n - length(out))))
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

# ---- the third reference: a plain-double difference of the rates -------------

# The whole right-hand side differenced in each ODE state entry, in plain double.
#
# This is a reference the forward tangent cannot be, and the reason is structural
# rather than a matter of accuracy. The environment holds its integrated state as
# a double store and takes it passively, so a tangent seeded on a soil moisture
# carries nothing: every soil column and every soil rate row of its Jacobian is
# exactly zero. A difference re-runs the forward balance, so it carries them.
#
# It is a legitimate referee here, and only here, because the soil balance has no
# supplied row in it. Where a row IS supplied -- the leaf's -- a difference of the
# step that consumes it returns zero whether the row is right, wrong or absent,
# which is why the block's rows are refereed by the tangent instead.
ladder_rhs_state_difference <- function(x, rel = 1e-6) {
  patch <- ladder_as_patch(x)
  n <- patch$ode_size
  state <- patch$ode_state
  time <- patch$ode_time
  at <- function(v) { patch$set_ode_state(v, time); patch$ode_rates }
  out <- matrix(0, n, n)
  for (j in seq_len(n)) {
    h <- max(abs(state[[j]]), 1) * rel
    up <- state; up[[j]] <- state[[j]] + h
    dn <- state; dn[[j]] <- state[[j]] - h
    out[, j] <- (at(up) - at(dn)) / (2 * h)
  }
  # Leave the fixture where it was found, and rebuild the auxiliaries a regime
  # report is read off.
  invisible(at(state))
  out
}

# The finer difference's own error, measured rather than chosen.
#
# `coarse` is taken at ten times `fine`'s step. A central difference is
# second-order while truncation dominates, so the coarse step's error is about a
# hundred times the fine one's and the gap between them is therefore the COARSE
# error. Dividing by that factor is what turns the gap into an estimate of the
# error of the step actually being used -- taking the gap itself as the tolerance
# leaves it loose by two orders, which would not catch a soil row wrong by ten.
ladder_difference_floor <- function(coarse, fine, order = 2) {
  scale <- max(max(abs(coarse)), max(abs(fine)), .Machine$double.xmin)
  gap <- max(abs(coarse - fine)) / scale
  max(gap / 10^order, 4 * .Machine$double.eps)
}

# The same parameter set with one species' strategy parameter moved.
#
# A round trip through the list is what writes it through: the accessors hand out
# copies, so the nested value has to be reassigned outward at each level.
ladder_perturbed_parameters <- function(species, index, name, value) {
  p <- ladder_parameters(species)
  strategies <- p$strategies
  pars <- strategies[[index]]$pars
  pars[[name]] <- value
  strategies[[index]]$pars <- pars
  p$strategies <- strategies
  p
}

# The whole right-hand side differenced in one strategy parameter by REBUILDING
# the strategy from its parameters, at the state the base fixture holds.
#
# This is the complement of the prepared-strategy difference, and the two are
# unbiased on disjoint sets. Perturbing a prepared strategy cannot reach a trait
# the leaf holds its own copy of, so that difference reads exactly zero there. A
# rebuild runs preparation, so the leaf receives the trait -- at the price of also
# moving the seed height, which the differentiated path imposes to zero.
#
# The price is not paid on the columns this is for. The seed height solves
# mass_live(h) = omega, which reads the allometric and tissue parameters and no
# leaf trait at all. So for exactly the columns the prepared perturbation is blind
# to, a rebuild carries no birth-size channel either, and each scheme is unbiased
# precisely where the other is not.
ladder_rate_difference_rebuilt <- function(base, index, name, rel = 1e-6) {
  state <- base$ode_state
  time <- base$ode_time
  value <- ladder_strategy_parameter(base, index, name)
  h <- max(abs(value), 1) * rel
  at <- function(v) {
    patch <- ladder_patch_two_by_two(
      cross = FALSE,
      parameters = ladder_perturbed_parameters(c("fast", "slow"), index, name, v))
    # The base fixture's own state, imposed rather than reconditioned, so nothing
    # but the parameter differs between the two evaluations.
    patch$set_ode_state(state, time)
    patch$ode_rates
  }
  (at(value + h) - at(value - h)) / (2 * h)
}

# The leaf's own traits, which no difference of the rates can referee.
#
# The leaf holds its own copy of each of these, taken when the strategy was
# prepared, and a rate evaluation does not push them back in. So perturbing the
# strategy's parameter leaves the physiology at the value it was prepared with and
# the difference comes back EXACTLY zero -- not small, zero -- whether the sweep's
# row is right, wrong or absent. The same shape as a supplied row one level down.
#
# Their referee is the leaf's own algebra at one solved operating point, which is
# rung 1 and is not built. Until it is, these columns are unrefereed by anything,
# and that is a wider gap than the soil channel the tangent misses.
ladder_leaf_own_traits <- function() {
  c("c", "b", "beta2", "g1_TF24", "a", "curv_fact_elec_trans",
    "curv_fact_colim", "root_c", "root_b")
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

# A rung blocked because the block cannot record at an active scalar is one red
# line, named once. Shared by every file that forms the block or consumes its
# rows, which is why it is here rather than in the first file that needed it.
ladder_block_or_skip <- function(patch, node = 1L) {
  out <- tryCatch(list(value = ladder_block_value_tf24(patch, node)),
                  error = function(e) e)
  if (inherits(out, "error")) {
    testthat::skip(paste("the cohort block does not record at an active scalar:",
                         conditionMessage(out)))
  }
  out
}

ladder_gradient_or_skip <- function(stand, ...) {
  result <- tryCatch(stand_gradient(stand, ...), error = function(e) e)
  if (inherits(result, "error")) {
    testthat::skip(paste("the sweep does not run on this stand:",
                         conditionMessage(result)))
  }
  result
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

# A parameter belongs on one of the lists below only when it has no live route at
# all; one whose birth-size channel alone is imposed to zero has a live column and
# belongs on the birth-size list further down.
#
# A declared zero is asserted and not merely tolerated, and each list is a claim
# about *why* its columns are zero that ladder_zero_cause() measures at the rate
# level: the cause names which rows of the right-hand side the parameter moves, so
# a column that is zero for a reason other than the one declared fails as loudly
# as an undeclared zero does.
#
# The lists are per level, because the levels have different zeros. A column zero
# in the cohort block can be live in the census and the reverse, so excusing one
# level's zero at the other is how an absent row gets excused.

# The dry bound of the leaf's feasible interval, and nothing inside it. The bound
# is the lesser of the stem's and the root's critical potentials, and at an
# interior operating point the solve never reaches it: moving one leaves profit,
# the collar and every per-layer draw bit-identical, so the row is exactly zero
# by complementary slackness rather than by omission.
#
# This is a property of the regime and not of the model. At a pinned operating
# point these two carry the whole row -- and where the root's potential wins the
# minimum, that row is exactly minus the unit vector in its own direction. The
# fixtures are pinned interior, so nothing here exercises it.
ladder_zero_at_an_interior_optimum <- function() {
  c("psi_crit", "root_psi_crit")
}

# The seed's mass and the accessory cost of a seed. They reach two rates and no
# others -- offspring production, and the survival-weighted offspring the census
# does not read -- and they reach both only through the sum (omega + a_f3), so
# the two are not separately identifiable even there. Neither accumulator is read
# by any equation or by any of the three metrics, so the column is exactly zero
# on this metric set and would be live on a fitness functional.
#
# omega has a second route, through birth size, and that one is the imposed zero
# below rather than this one.
ladder_zero_outside_the_metric_support <- function() {
  c("omega", "a_f3")
}

# The rates a parameter outside the metric support is allowed to move. Naming
# them is what makes the class a measurement: a third rate would mean the
# census's silence about the column is wrong.
ladder_reproductive_rates <- function() {
  c("fecundity", "offspring_produced_survival_weighted")
}

# Reach the census through the introduction boundary and through nothing else.
# The block is one individual's physiology and reads none of them: the
# establishment probability, its decay in patch age and the initial reserve are
# evaluated once per introduction, at the seed size, in the field the newcomer is
# placed in. So an exact zero in the BLOCK is correct here and an exact zero in
# the census is not -- which is the whole content of the class, and is why the
# two levels carry different declared lists.
#
# The other half of the claim is asserted where it can be: the floor's census
# classification puts all three in the live class, and rung 5's single-channel
# probes take a_st3 and recruitment_decay by name.
ladder_zero_outside_the_cohort_block <- function() {
  c("a_d0", "a_st3", "recruitment_decay")
}

ladder_zero_by_construction <- function() {
  c(ladder_zero_at_an_interior_optimum(),
    ladder_zero_outside_the_metric_support())
}

# Which cause a declared zero is claimed to have. Named so a check can measure
# the claim rather than only record it.
ladder_zero_cause <- function(name) {
  bare <- ladder_bare_traits(name)
  if (bare %in% ladder_zero_at_an_interior_optimum()) {
    "interior optimum: moves no rate at all"
  } else if (bare %in% ladder_zero_outside_the_metric_support()) {
    "outside the metric support: moves the two offspring rates"
  } else {
    NA_character_
  }
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

# The trajectory reference: a tangent run of the same solve, stepped at the sizes
# the run recorded. `direction` carries one weight per gradient column, so a unit
# vector asks for one column and a mixed one for a contraction.
#
# It returns the metrics the replay itself reached as well as the derivative,
# because a reference whose value disagrees with the model is a reference to a
# different trajectory and every derivative taken from it is that trajectory's.
ladder_trajectory_tangent <- function(stand, direction) {
  ladder_trajectory_tangent_tf24(stand, direction)
}

# One column's disagreement with a reference, per metric row and normalised per
# row.
#
# Normalising a whole column by its largest row is what a single number invites and
# it hides the small rows: leaf area and above-ground mass are of order one while
# stem area is of order 1e-4, so a stem-area row wrong by three per cent reads as
# 1e-6 of the column's peak. Report 08 §4 asks for residuals per block for exactly
# this reason, and a column is a block here.
#
# `trajectory` is the column with the census's direct term removed. Where the total
# is a near-cancellation of the two terms -- which it is for the allometric
# constants -- the relative error against the total is an amplified view of the
# error in the trajectory term, and the second number is the unamplified one.
ladder_column_residual <- function(got, reference, trajectory = NULL) {
  scale <- pmax(abs(reference), max(abs(reference)) * 1e-12)
  out <- list(per_row = max(abs(got - reference) / scale))
  if (!is.null(trajectory)) {
    out$per_row_trajectory <-
      max(abs(got - reference) / pmax(abs(trajectory),
                                      max(abs(trajectory)) * 1e-12))
    out$amplification <- max(abs(trajectory) / scale)
  }
  out
}

# A unit direction in trait space, named as the gradient's columns are.
ladder_trait_direction <- function(columns, name) {
  at <- match(name, columns)
  stopifnot(!is.na(at))
  replace(numeric(length(columns)), at, 1)
}

# How far the replay's own census lands from the run's. This is the trajectory
# checks' floor: the reference replays step sizes rather than times, so its state
# at the end is the run's state only up to the arithmetic of re-stepping, and
# every derivative it reports is a derivative at the state it actually reached.
ladder_replay_floor <- function(stand, value) {
  reached <- ladder_trajectory_tangent(
    stand, numeric(length(census_trait_names_tf24(stand))))$value
  max(abs(reached - value) / abs(value))
}

# What the trajectory sweep and the trajectory tangent agree to, for a column
# whose route to the census does not run through the newcomer's own leaf area.
# Measured worst 6.0e-05 over all 88 columns of a stand carrying one cohort per
# species, and 3.4e-05 over the four-node stand's columns outside the class below.
#
# The two parameters of the inflow condition itself are held here, and were not
# before: recruitment_decay sat at 1.9e-02 and a_d0 at 2.0e-03 while the boundary
# node's density adjoint was accumulated and never transposed through the condition
# that sets it. They now measure 2.9e-12 and 2.1e-05.
ladder_trajectory_agreement <- function() 3e-04

# And what they agree to for the two allometric constants once a species carries
# more than one cohort. Measured worst 3.4e-03, and it is an open defect rather
# than a floor.
#
# What is known about it, because each step was measured rather than argued:
#
#   It is per species and it switches on with that species' SECOND cohort. With
#   one cohort each, a_l2 measures 6.0e-05 and 1.4e-06; introduce a second cohort
#   of species 1 only and species 1 jumps to 3.4e-03 while species 2 stays at
#   1.6e-06; introduce a second of species 2 and species 2 jumps to 2.8e-03.
#   Further cohorts do not add to it -- 3.4e-03, 1.8e-03, 2.2e-03 -- so it appears
#   once and saturates.
#
#   It is not the inflow condition. Closing that transpose moved recruitment_decay
#   by ten orders and a_d0 by two, and left this unchanged at 3.05e-03 against
#   3.06e-03.
#
#   It is not the state a rebuild linearises at. Loading a recorded state without
#   re-evaluating the condition in the field the nodes were rated in was a second
#   defect, worth 1.2 per cent in the boundary density by the end of a run; fixing
#   it improved these columns tenfold on a one-cohort stand and did not move them
#   here.
#
#   It is not k_I, which shares the field reduction and measures 3.3e-07.
#
#   It is not the census's own reading of the traits. That term is now reported on
#   its own and refereed against a plain-double difference of the census: they
#   agree to 5.1e-09 over the matrix, which is the difference's round-off. So the
#   whole of the disagreement is in the trajectory term.
#
#   It is not the reduction's own parameter rows. Rung 3 forms them entry by entry
#   on a two-cohort fixture, against a forward Jacobian, and takes a_l1 and a_l2 by
#   name among the reduction-borne columns.
#
#   It is not the rows an introduction carries through untouched. Contracting the
#   whole widened state rather than copying those rows and contracting the
#   newcomer's leaves every figure here bit-identical, so they are the exact
#   identity the copy assumed.
#
#   It is not the introduction boundary at all, which is the result that cost the
#   most to get and is worth stating plainly. Three checks now cover it and all
#   three are clean: the introduction map's whole Jacobian agrees forward against
#   reverse to 7.5e-16 over every cell at both widenings, the rows an introduction
#   carries through are bit-identical to the recorded state, and the rows it writes
#   are bit-identical to the boundary node the patch held. The newcomer's own
#   sensitivity to these two traits is 4.0e-07 and 1.9e-05, because the seed height
#   is imposed passive and only the seed's physiology is left.
#
# So what remains is a species carrying two cohorts rather than one, with the
# introduction that necessarily accompanies it excluded. Rung 3 forms the per-stage
# Jacobian entry by entry on a two-cohort patch and finds it correct, so the
# remaining difference between the two is the TIME dimension.
#
# It has a per-step character and is not a clean lost term. Tightening the ODE
# tolerance from 1e-04 to 1e-10 takes the run from 98 steps to 635 and the residual
# from 3.3e-02 to 1.3e-02 -- six and a half times the steps for two and a half
# times the residual, near h^0.6. The sweep and the tangent replay the SAME
# recorded step sizes, so a truncation error would cancel between them exactly;
# what does not cancel is a linearisation taken at a point the two do not share.
# k_I moves over the same range from 8.8e-08 to 9.3e-09, so the effect is not
# uniform across traits either.
#
#   It is not the stage recursion, and that one was settled by injection rather
#   than by reading. Replacing the sweep's scatter over every earlier stage with
#   the immediate predecessor alone -- report 08 §7's named fault, built into
#   odelia and measured against this pair -- moves the ONE-cohort control from
#   2.2e-07 to 7.0e-03, a margin of thirty-one thousand, and saturates the longer
#   fixtures at 1.0. So the check is live and enormously sensitive to a lost stage
#   term, and the real defect is not one: a lost term raises the one-cohort control
#   to the same order as the two-cohort case, while the defect leaves it at 2.2e-07.
#   The recursion is correct by inspection too -- it scatters h * a_im over every
#   m < i on the dense Cash-Karp rows, accumulates the state channel at every stage,
#   and seeds the rate accumulators from the output weights, leaving the two whose
#   weights vanish empty.
#
#   It is not the size-distribution quadrature's interior interval either, which is
#   the first thing a second cohort creates that a single one does not. Moving the
#   second introduction from 0.02 to 0.38 of a four-tenths run widens that interval
#   from 0.02 to 0.38 and the residual FALLS, 7.6e-03 to 4.7e-04. It does not track
#   the interval.
#
#   It is not the stage state and aux restore, which is the one per-step path rung 3
#   does not exercise. Instrumenting that boundary -- recomputing the rates after
#   the restore and comparing against the stage the forward pass computed -- gives a
#   worst relative difference of ZERO over every stage of every step, on one cohort
#   and on two. The restore is faithful.
#
#   It is not the light field's discretisation, and this one is worth its own line
#   because report 03 §3.3 asks for the measurement and nobody had taken it. Built
#   at 33, 65 and 129 knots, a_l2's residual is 1.4371e-03 at all three -- identical
#   to five figures, while the forward model itself moves. k_I's, by contrast, falls
#   from 2.22e-10 to 4.58e-11 over the same range. So the passive knot-position
#   treatment DOES converge with knot density for a field-borne trait, which is
#   report 03's own falsifier answered in its favour, and a_l2's residual is not a
#   field quantity at all.
#
# What it tracks is the time the stand spends carrying two cohorts, linearly and at
# about 0.02 per unit time: 7.6e-03 over 0.38, 4.0e-03 over 0.20, 4.7e-04 over 0.02.
# Three cohorts land where two of the same duration would, so it turns on at the
# second cohort and does not scale with the count. So it is accumulated per step
# while a species carries more than one cohort -- which is what makes the remaining
# suspects the per-step work that rung 3 does NOT exercise, the stage state and aux
# restore in particular: step_adjoint reloads each stage with
# set_ode_state_and_field and then set_ode_aux from the values captured on the
# forward pass, and both vectors are wider with two cohorts than with one. Rung 3
# takes its Jacobian on a patch whose aux are current, so nothing in this ladder
# covers that restore.
#
# Two things about how it is measured, both of which were wrong before and are
# worth keeping right. A column is normalised PER ROW, because leaf area and
# above-ground mass are of order one while stem area is of order 1e-4 and a
# stem-area row wrong by three per cent reads as 1e-06 of the column's peak. And
# where the total is a near-cancellation of the direct and trajectory terms -- on a
# stand whose species carries two cohorts the ratio reaches -16.8 on stem area --
# the relative error against the total is an amplified view of the trajectory
# term's, so ladder_column_residual() reports both and the amplification beside
# them.
#
# On this coordinate a second cohort cannot exist without an introduction after
# t = 0, because tied birth dates are refused, so "per introduction" and "per
# interior trapezium of the reduction" are not separable by scheduling alone.
# Seeding the cohorts instead does not separate them either: a run resumed from a
# populated state lands far outside this ladder's declared regime -- relative
# reserve 0.44 to 0.73 against the band 0.02 to 0.30, gate slope 0.010 against the
# floor 0.4 -- so such a run is invalid as an instrument rather than failing, and
# its own disagreements are larger again. The check that would separate them is the
# introduction map's own Jacobian, formed entry by entry against a tangent, which
# is the one unit-level object rung 5 does not yet have.
ladder_introduction_residual <- function() 4e-03

# Which shortlisted parameters are held to the looser bound above: the two
# constants that set leaf area from height, and no others.
ladder_introduction_borne_traits <- function() c("a_l1", "a_l2")

# What the transpose and the two references agree to on the soil rows, which is
# the tightest the leaf's supplied rows allow: a difference of the block re-solves
# the collar and the transpose grafts, so the gap is the rank-two factorisation's
# own fit step rather than a missing channel.
#
# It was 1.06e-05 while the uptake reduction's transpose scattered onto the
# introduced nodes only, leaving the boundary node's own draw -- the distribution's
# lower grid point -- reaching no soil row. Both reductions now carry that row and
# the residual is 5.38e-08 on the four-node patch, 5.07e-08 crossed and 4.67e-08
# on the one-cohort patch.
ladder_soil_row_agreement <- function() 8e-08

# A parameter the boundary cannot answer for. A refusal, never a number, and an
# unknown parameter must refuse by name rather than return anything at all.
ladder_refused_by_name <- function() {
  c("p_50")
}

# The band a registered parameter that reaches no equation comes back in. It is
# round-off rather than exact zero, so the check is two-sided and not a test
# against zero.
ladder_roundoff_band <- function() c(1e-22, 1e-16)
