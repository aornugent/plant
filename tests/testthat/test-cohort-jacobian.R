# An independent finite-difference reference for one TF24 cohort's rate vector
# at one solver stage, for a future reverse-mode sensitivity pass to be checked
# against. See R/cohort_jacobian.R.
#
# The claims here are of three kinds, in increasing strength:
#
#  1. the harness isolates one cohort -- a perturbation moves what it should and
#     nothing else, and leaves no residue;
#  2. the difference quotient has a step-size plateau, so it is a derivative
#     rather than truncation error or cancellation;
#  3. the sparsity pattern is exact. Some state components cannot influence some
#     rates at all, and those entries are asserted identically zero rather than
#     small. That is a statement about the model, and its transpose is what an
#     adjoint has to satisfy.

# One cohort part-way through a stand's life, in a real canopy and a real soil
# column, with the environment frozen at that instant. Built once per file.
cohort_jacobian_fixture <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    env <- Environment("TF24")
    env$set_soil_number_of_depths(5)
    env$set_soil_water_state(rep(0.3, 5))

    p <- scm_base_parameters("TF24")
    p$max_patch_lifetime <- 2
    p <- add_strategies(p, trait_matrix(0.0825, "lma"))

    scm <- SCM("TF24", "TF24_Env")(p, env, control())
    scm$run()

    patch <- scm$patch
    species <- patch$species[[1]]
    # A node in the interior of the size distribution: neither the boundary
    # node nor the tallest, so no state sits on a clamp.
    node <- species$node_at(round(species$size / 2))
    environment <- patch$environment
    pr_patch_survival <- patch$pr_survival(patch$time)
    # Shared by the two step-size tests; ~450 rate evaluations, about 1.5 s.
    sweep <- cohort_jacobian_sweep(node, environment, pr_patch_survival,
                                   steps = 10^seq(-11, -2, by = 0.25))
    cached <<- list(node = node, environment = environment,
                    pr_patch_survival = pr_patch_survival, sweep = sweep)
    cached
  }
})

test_that("one cohort can be isolated and perturbed without residue", {
  f <- cohort_jacobian_fixture()
  node <- f$node

  x0 <- node$ode_state
  r0 <- cohort_rate_vector(node, f$environment, f$pr_patch_survival)

  # Repeating the evaluation reproduces it bit for bit: the leaf sub-model is
  # held on the strategy and reused across calls, so this is not free.
  expect_identical(cohort_rate_vector(node, f$environment, f$pr_patch_survival),
                   r0)

  # And the evaluation is path-independent: going somewhere else in state space
  # and coming back gives the same numbers, bit for bit.
  moved <- x0
  moved[[1]] <- moved[[1]] * 1.05
  node$ode_state <- moved
  cohort_rate_vector(node, f$environment, f$pr_patch_survival)
  node$ode_state <- x0
  expect_identical(cohort_rate_vector(node, f$environment, f$pr_patch_survival),
                   r0)

  # Forming the whole Jacobian leaves the node where it was found.
  cohort_rate_jacobian(node, f$environment, f$pr_patch_survival)
  expect_identical(node$ode_state, x0)
  expect_identical(node$ode_rates, unname(r0[node$ode_names]))
})

test_that("states nothing reads back cannot move any rate at all", {
  # The strongest form of the sparsity claim, and the one that does not depend
  # on the finite-difference step: move these states by an amount that is not
  # small, and the entire rate vector has to come back bit-identical. A
  # tolerance would not distinguish "does not depend on" from "depends weakly".
  f <- cohort_jacobian_fixture()
  node <- f$node
  x0 <- node$ode_state
  on.exit(node$ode_state <- x0, add = TRUE)

  r0 <- cohort_rate_vector(node, f$environment, f$pr_patch_survival)
  names(x0) <- node$ode_names

  for (nm in c("fecundity", "area_heartwood", "mass_heartwood")) {
    x <- x0
    x[[nm]] <- x0[[nm]] + 1
    node$ode_state <- unname(x)
    expect_identical(
      cohort_rate_vector(node, f$environment, f$pr_patch_survival), r0,
      info = nm)
  }

  # Mortality is different: it IS read, but only through a finiteness test in
  # mortality_dt and through exp(-mortality) in the node's offspring rate. So
  # every rate but that one is invariant to it, exactly.
  x <- x0
  x[["mortality"]] <- x0[["mortality"]] + 0.5
  node$ode_state <- unname(x)
  r <- cohort_rate_vector(node, f$environment, f$pr_patch_survival)
  moved <- "offspring_produced_survival_weighted"
  expect_identical(r[names(r) != moved], r0[names(r0) != moved])
  expect_equal(r[[moved]], r0[[moved]] * exp(-0.5))
})

test_that("the cohort Jacobian's structural zeros are exactly zero", {
  f <- cohort_jacobian_fixture()
  pattern <- cohort_jacobian_pattern(f$node)

  # 58 of the 78 entries are structurally zero at five soil layers.
  expect_equal(sum(pattern), 58L)

  jacobian <- cohort_rate_jacobian(f$node, f$environment, f$pr_patch_survival)
  expect_identical(dim(jacobian), dim(pattern))
  expect_identical(jacobian[pattern], rep(0, sum(pattern)))

  # And the entries the pattern claims are live are in fact live, so the
  # pattern is not passing by being over-broad.
  expect_true(all(jacobian[!pattern] != 0))
})

test_that("the cohort Jacobian reproduces the derivatives known in closed form", {
  # Four entries can be differentiated by hand from the model equations, which
  # is the only check available that does not go through the same finite
  # difference it is checking. Two are in the height column, one in the
  # mortality column and one in the storage column, so all three live columns
  # are covered.
  f <- cohort_jacobian_fixture()
  node <- f$node
  pars <- node$individual$strategy$pars
  state <- stats::setNames(node$ode_state, node$ode_names)
  rates <- cohort_rate_vector(node, f$environment, f$pr_patch_survival)

  h <- state[["height"]]
  eta_c <- 1 - 2 / (1 + pars$eta) + 1 / (1 + 2 * pars$eta)
  area_leaf <- (h / pars$a_l1)^(1 / pars$a_l2)
  darea_leaf_dh <- area_leaf / (pars$a_l2 * h)

  jacobian <- cohort_rate_jacobian(node, f$environment, f$pr_patch_survival)

  # area_heartwood_dt = k_s * theta * area_leaf(h)
  expect_equal(jacobian["area_heartwood", "height"],
               pars$k_s * pars$theta * darea_leaf_dh, tolerance = 1e-6)

  # mass_heartwood_dt = k_s * theta * area_leaf(h) * h * eta_c * rho
  expect_equal(jacobian["mass_heartwood", "height"],
               pars$k_s * pars$theta * eta_c * pars$rho *
                 (darea_leaf_dh * h + area_leaf),
               tolerance = 1e-6)

  # offspring_dt is proportional to exp(-mortality).
  expect_equal(jacobian["offspring_produced_survival_weighted", "mortality"],
               -rates[["offspring_produced_survival_weighted"]],
               tolerance = 1e-6)

  # mortality_dt = d_I + a_dG1 * exp(-a_dG2 * S / S_max),
  # S_max = a_st1 * theta * area_leaf(h) * h * eta_c * rho.
  storage_max <- pars$a_st1 * pars$theta * area_leaf * h * eta_c * pars$rho
  relative_reserves <- state[["storage"]] / storage_max
  expect_lt(relative_reserves, 1)  # off the r = 1 clamp, or the entry is 0
  expect_equal(jacobian["mortality", "storage"],
               -pars$a_dG1 * pars$a_dG2 *
                 exp(-pars$a_dG2 * relative_reserves) / storage_max,
               tolerance = 1e-6)
})

test_that("the finite-difference step has a plateau with both edges", {
  # Do not pick a step. The difference quotient is only a derivative over a
  # range of steps, bounded below by cancellation against the leaf sub-model's
  # own solver tolerances and above by curvature; outside it the method returns
  # something plausible and wrong. Establish the range, then take its centre.
  f <- cohort_jacobian_fixture()
  sweep <- f$sweep
  plateau <- cohort_jacobian_plateau(sweep, tol = 1e-6)
  physiology <- plateau[plateau$rate != "log_density", ]

  # Every live entry that depends on the physiology has one, and every one of
  # them has a lower edge interior to the swept range.
  live <- !cohort_jacobian_pattern(f$node)
  expect_equal(nrow(physiology), sum(live[rownames(live) != "log_density", ]))
  expect_true(all(physiology$lower > 1e-11))

  # Individually, an entry flat enough over the whole grid need not show an
  # upper edge -- the offspring rate is close to linear in the mortality state,
  # so its truncation error stays below the tolerance out to 1e-2. The claim
  # that matters is about the range they share: it is bounded on both sides,
  # wide, and contains the step the harness defaults to.
  common <- c(max(physiology$lower), min(physiology$upper))
  expect_gt(common[[1]], 1e-11)
  expect_lt(common[[2]], 1e-2)
  expect_gt(log10(common[[2]] / common[[1]]), 1.5)
  expect_lte(common[[1]], formals(cohort_rate_jacobian)$step)
  expect_gte(common[[2]], formals(cohort_rate_jacobian)$step)
})

test_that("the density rate does not share the plateau, and why", {
  # log_density_dt is built from a backward difference of the height growth
  # rate at a fixed absolute step (control()$node_gradient_eps = 1e-6), so
  # differencing it again inherits that difference's roundoff. Measured on this
  # fixture, the density rate's relative noise floor is 4e-8 against 5e-10 to
  # 4e-9 for the rates that come straight out of the physiology.
  #
  # This is the failure the sweep exists to catch: at a step chosen for the
  # other rates, these two entries are not converged, and they do not announce
  # it -- they return a smooth, plausible number.
  f <- cohort_jacobian_fixture()
  sweep <- f$sweep
  default_step <- formals(cohort_rate_jacobian)$step

  fine <- cohort_jacobian_plateau(sweep, tol = 1e-6)
  fine <- fine[fine$rate == "log_density", ]
  expect_false(any(fine$lower <= default_step & fine$upper >= default_step))

  # A hundred times looser, and the density rate does have a plateau -- above
  # the one everything else uses, on both entries.
  coarse <- cohort_jacobian_plateau(sweep, tol = 1e-4)
  coarse <- coarse[coarse$rate == "log_density", ]
  expect_equal(nrow(coarse), 2L)
  expect_true(all(coarse$lower > default_step))
})
