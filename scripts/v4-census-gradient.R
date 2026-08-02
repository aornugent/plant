# V4: TF24 census and R0 at max_patch_lifetime = 105.32, resident, against a
# re-run finite difference on the identical resolved schedule.
#
# Run with R_LIBS_USER pointing at the odelia install and the plant tree loaded
# from source:
#   R_LIBS_USER=/home/user/lib-p3-int Rscript scripts/v4-census-gradient.R
#
# The harness prints which traits it covered and which it did not. A central
# difference over all 44 of TF24's differentiable parameters is 88 production
# runs; the subset below is nine traits, so eighteen, chosen so each discriminates
# a different channel, and the traits left out are named rather than silently
# dropped.

suppressMessages(library(odelia))
pkgload::load_all(".", quiet = TRUE)

MAX_PATCH_LIFETIME <- 105.32

# The establishment gate takes its closed arm on 23.1% of boundary-node stage
# evaluations, confined to this window. A trait perturbation of ~1e-9 flips it,
# so a re-run finite difference inside the window disagrees with the adjoint for
# a reason that is neither's error. V4 therefore differences a census evaluated
# outside it.
ESTABLISHMENT_WINDOW <- c(3.222267, 8.544184)

# At the end of a run the bottom of the size distribution has underflowed to
# survival exactly zero, so the closing trapezium down to the boundary node
# contributes exactly nothing and a census differenced at T carries no
# boundary-node channel at all. The harness reports the live/dead split so that
# reads as a property of the state rather than as a passing check.

# Each trait is here because it reaches the census through a channel no other
# trait in the list reaches it through.
V4_TRAITS <- list(
  lma = "leaf mass per area: the leaf construction cost and the leaf-area
         allometry's mass side. Moves every metric through mass_leaf directly and
         through the carbon balance, so a zero row here means nothing is wired.",
  rho = "wood density: enters mass_sapwood and mass_bark only through the
         allometry, and the carbon balance through sapwood cost. The one trait
         that separates the mass metric from the area metrics.",
  hmat = "height at maturation: moves the allocation to reproduction, so it is
          the trait that reaches R0 strongly and the standing metrics weakly.
          Discriminates the two functionals' different sensitivity structures.",
  theta = "sapwood area per leaf area: reaches area_stem through area_sapwood and
           area_bark and reaches nothing else in the allometry. The cleanest test
           that the basal-area metric is not just tracking leaf area.",
  a_l1 = "the height-leaf area allometry's scale: moves every cohort's height at
          a given leaf area, so it moves the quadrature grid. The trait that
          exposes a dropped trapezium-weight term.",
  k_I = "light extinction: reaches the census only through the canopy, so it is
         the resident channel (step (c)) and is exactly zero on an invasion
         gradient. Distinguishes the two passes.",
  a_dG1 = "growth-dependent mortality: reaches the census only through survival
           and hence through log density, with no allometric route at all.",
  K_s = "stem-specific conductivity: reaches the census only through the leaf
         solve's operating point, so it is the one trait in this list whose row
         is zero if Leaf::input_adjoints is not connected.",
  psi_crit = "the critical water potential: the hydraulic cost term, and the
              trait most likely to move the establishment gate. Included so the
              window avoidance is exercised rather than assumed."
)

# Step sizes. A central difference of a trait resolves the census to about four
# digits at h ~ 1e-5 of the trait's own scale; 1e-9 is the scale at which the
# establishment gate flips, so nothing here goes near it.
FD_RELATIVE_STEP <- 1e-5

build_stand <- function(lifetime = MAX_PATCH_LIFETIME) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(0.0825, "lma"))
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

v4 <- function() {
  p <- build_stand()
  base <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  base$refine_schedule()
  schedule <- base$node_schedule$all_times
  ode_times <- base$ode_times

  census_time <- base$time
  in_window <- census_time >= ESTABLISHMENT_WINDOW[1] &&
    census_time <= ESTABLISHMENT_WINDOW[2]
  message(sprintf(
    "census evaluated at t = %.4f; establishment window is [%.4f, %.4f]: %s",
    census_time, ESTABLISHMENT_WINDOW[1], ESTABLISHMENT_WINDOW[2],
    if (in_window) "INSIDE - the finite difference is not the reference here"
    else "outside — the finite difference is a valid reference"))
  if (in_window) {
    stop("V4's census time is inside the establishment window; move it out ",
         "rather than mollifying the model")
  }

  sp <- base$patch$species[[1]]
  nm <- sp$new_node$ode_names
  state <- matrix(c(sp$ode_state, sp$new_node$ode_state), nrow = length(nm),
                  dimnames = list(nm, NULL))
  density <- exp(state["log_density", ])
  message(sprintf(
    "%d of %d cohorts at density exactly zero; boundary node density %g: %s",
    sum(density == 0), length(density), density[[length(density)]],
    if (density[[length(density)]] > 0 && density[[length(density) - 1]] > 0)
      "the closing trapezium is live"
    else "the closing trapezium contributes nothing at this state"))

  adjoint <- stand_gradient(base, traits = names(V4_TRAITS))

  all_traits <- census_trait_names_tf24(base)
  covered <- names(V4_TRAITS)
  message("traits covered by V4: ", paste(covered, collapse = ", "))
  message("traits NOT covered (", length(setdiff(all_traits, covered)), " of ",
          length(all_traits), "): ",
          paste(setdiff(all_traits, covered), collapse = ", "))
  for (tr in covered) {
    message("  ", tr, ": ", gsub("\\s+", " ", V4_TRAITS[[tr]]))
  }

  fd <- matrix(NA_real_, nrow = length(adjoint$value), ncol = length(covered),
               dimnames = list(names(adjoint$value), covered))
  r0_fd <- stats::setNames(rep(NA_real_, length(covered)), covered)
  for (tr in covered) {
    base_value <- p$strategies[[1]]$pars[[tr]]
    h <- FD_RELATIVE_STEP * max(abs(base_value), 1)
    out <- lapply(c(-1, 1), function(sign) {
      pp <- p
      pp$strategies[[1]]$pars[[tr]] <- base_value + sign * h
      scm <- run_on_schedule(pp, schedule, ode_times)
      list(census = stand_census(scm),
           r0 = scm$net_reproduction_ratios[[1]])
    })
    fd[, tr] <- (out[[2]]$census - out[[1]]$census) / (2 * h)
    r0_fd[[tr]] <- (out[[2]]$r0 - out[[1]]$r0) / (2 * h)
  }

  list(value = adjoint$value,
       adjoint = adjoint$gradient,
       finite_difference = fd,
       r0_finite_difference = r0_fd,
       control = adjoint$control,
       traits_covered = covered,
       traits_not_covered = setdiff(all_traits, covered),
       census_time = census_time)
}

if (identical(environment(), globalenv())) {
  res <- v4()
  print(res$adjoint)
  print(res$finite_difference)
  rel <- abs(res$adjoint - res$finite_difference) /
    pmax(abs(res$finite_difference), 1e-30)
  message("worst relative disagreement: ", format(max(rel)))
  message("peak memory (MB): ",
          format(sum(gc()[, "max used"] * c(8, 56)) / 2^20))
}
