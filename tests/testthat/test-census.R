# The census reduction and the R-facing gradient entry point.

# TF24's allometry, written out from its equations rather than called through
# the strategy, so the value gate does not check the reduction against itself.
tf24_allometry_r <- function(pars, eta) {
  eta_c <- 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta)
  function(height, area_heartwood, mass_heartwood) {
    area_leaf <- (height / pars$a_l1)^(1 / pars$a_l2)
    area_sapwood <- area_leaf * pars$theta
    area_bark <- pars$a_b1 * area_leaf * pars$theta
    mass_leaf <- area_leaf * pars$lma
    mass_sapwood <- area_sapwood * height * eta_c * pars$rho
    mass_bark <- area_bark * height * eta_c * pars$rho
    list(leaf_area = area_leaf,
         mass_above_ground = mass_leaf + mass_bark + mass_sapwood +
           mass_heartwood,
         area_stem = area_bark + area_sapwood + area_heartwood)
  }
}

# Every cohort's state, boundary node last, from the ODE state the species and
# its boundary node write.
species_state_r <- function(species) {
  names <- species$new_node$ode_names
  stride <- length(names)
  flat <- c(species$ode_state, species$new_node$ode_state)
  s <- matrix(flat, nrow = stride, dimnames = list(names, NULL))
  data.frame(height = s["height", ],
             area_heartwood = s["area_heartwood", ],
             mass_heartwood = s["mass_heartwood", ],
             log_density = s["log_density", ])
}

trapezium_r <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1))) / 2
}

# The census in R: the grid runs from the boundary node upwards, so the state
# columns are reversed and the boundary node (last) becomes the first grid point.
census_r <- function(species, pars, eta, include_boundary = TRUE) {
  st <- species_state_r(species)
  n <- nrow(st)
  order <- if (include_boundary) c(n, seq(n - 1, 1)) else seq(n - 1, 1)
  st <- st[order, , drop = FALSE]
  psi <- tf24_allometry_r(pars, eta)(st$height, st$area_heartwood,
                                     st$mass_heartwood)
  density <- exp(st$log_density)
  vapply(psi, function(p) trapezium_r(st$height, density * p), numeric(1))
}

solved_stand <- function(lifetime = 20) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), Control())
  scm$run()
  scm
}

test_that("G1: the census value matches an independent R reduction", {
  scm <- solved_stand()
  species <- scm$patch$species[[1]]
  strategy <- scm$parameters$strategies[[1]]
  expected <- census_r(species, strategy$pars, strategy$pars$eta)
  got <- stand_census(scm)
  expect_equal(names(got), c("leaf_area", "mass_above_ground", "area_stem"))
  expect_equal(unname(got), unname(expected[names(got)]), tolerance = 1e-12)
})

test_that("G2: the boundary node is in the reduction", {
  # At a lifetime where the bottom of the size distribution is alive. Where the
  # boundary node's own density has underflowed to zero, so has its neighbour's,
  # the closing trapezium contributes exactly nothing and this comparison is
  # vacuous -- which is the case at t = 5 and t = 20.
  scm <- solved_stand(12)
  species <- scm$patch$species[[1]]
  strategy <- scm$parameters$strategies[[1]]
  st <- species_state_r(species)
  # The premise: the closing interval has live density at both ends.
  expect_true(exp(st$log_density[nrow(st)]) > 0)
  expect_true(exp(st$log_density[nrow(st) - 1]) > 0)
  with_boundary <- census_r(species, strategy$pars, strategy$pars$eta, TRUE)
  without <- census_r(species, strategy$pars, strategy$pars$eta, FALSE)
  # A reduction that starts at the smallest cohort drops the interval down to
  # the boundary node.
  expect_true(all(abs(with_boundary - without) / abs(with_boundary) > 1e-6))
  expect_equal(unname(stand_census(scm)), unname(with_boundary),
               tolerance = 1e-12)
})

test_that("G3: the census moves through the trapezium weights", {
  scm <- solved_stand(12)
  species <- scm$patch$species[[1]]
  strategy <- scm$parameters$strategies[[1]]
  st <- species_state_r(species)
  psi_at <- tf24_allometry_r(strategy$pars, strategy$pars$eta)

  # Perturb one interior cohort's height and difference the census. The
  # integrand-only census holds the grid fixed at the unperturbed heights, so the
  # difference between the two differences is the quadrature-weight term.
  k <- floor(nrow(st) / 4)
  eps <- 1e-6
  full <- function(h) {
    hh <- st$height
    hh[k] <- h
    order <- c(nrow(st), seq(nrow(st) - 1, 1))
    p <- psi_at(hh[order], st$area_heartwood[order], st$mass_heartwood[order])
    trapezium_r(hh[order], exp(st$log_density[order]) * p$leaf_area)
  }
  integrand_only <- function(h) {
    hh <- st$height
    hh[k] <- h
    order <- c(nrow(st), seq(nrow(st) - 1, 1))
    p <- psi_at(hh[order], st$area_heartwood[order], st$mass_heartwood[order])
    trapezium_r(st$height[order], exp(st$log_density[order]) * p$leaf_area)
  }
  d_full <- (full(st$height[k] + eps) - full(st$height[k] - eps)) / (2 * eps)
  d_part <- (integrand_only(st$height[k] + eps) -
               integrand_only(st$height[k] - eps)) / (2 * eps)
  # The weight term is the whole of this derivative and then some: it is larger
  # than the total and of the opposite sign to the integrand's own contribution.
  expect_true(abs(d_full - d_part) > abs(d_full))
  expect_true(sign(d_full) != sign(d_part))

  # The seed the reverse pass is given carries the whole derivative, weights
  # included, so it agrees with the full difference and not with the partial one.
  seed <- stand_census_state_adjoint(scm)
  stride <- length(scm$patch$species[[1]]$new_node$ode_names)
  col <- (k - 1) * stride + 1
  expect_equal(unname(seed["leaf_area", col]), d_full, tolerance = 1e-4)
})

test_that("G5: the entry point refuses to compare across two Controls", {
  a <- list(gradient = matrix(1, 1, 1, dimnames = list("leaf_area", "lma")),
            control = c(GSS_tol_abs = 1e-1, ci_abs_tol = 1e-3,
                        node_gradient_eps = 1e-6, schedule_eps = 1e-3))
  b <- a
  b$control[["schedule_eps"]] <- 1e-4
  expect_error(stand_gradient_compare(a, b), "schedule_eps")
  expect_equal(stand_gradient_compare(a, a),
               matrix(0, 1, 1, dimnames = list("leaf_area", "lma")))
})

test_that("G6: no census metric has an all-zero state sensitivity", {
  scm <- solved_stand()
  seed <- stand_census_state_adjoint(scm)
  expect_equal(nrow(seed), 3L)
  # Every metric is built from height, and every cohort's log density multiplies
  # it, so both state families must move all three metrics.
  stride <- length(scm$patch$species[[1]]$new_node$ode_names)
  n_node <- scm$patch$species[[1]]$size
  height_cols <- (seq_len(n_node) - 1) * stride + 1
  density_cols <- (seq_len(n_node) - 1) * stride + stride
  for (m in rownames(seed)) {
    expect_true(any(seed[m, height_cols] != 0))
    expect_true(any(seed[m, density_cols] != 0))
  }
})

test_that("the Control a gradient is taken at is the four that move it", {
  scm <- solved_stand()
  expect_equal(names(gradient_control(scm)),
               c("GSS_tol_abs", "ci_abs_tol", "node_gradient_eps",
                 "schedule_eps"))
  ctrl <- Control()
  expect_equal(unname(gradient_control(scm)),
               c(ctrl$GSS_tol_abs, ctrl$ci_abs_tol, ctrl$node_gradient_eps,
                 ctrl$schedule_eps))
})
