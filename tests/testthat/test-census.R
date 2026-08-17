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
# its boundary node write. The birth date is not ODE state and is read off the
# nodes themselves.
species_state_r <- function(species) {
  names <- species$new_node$ode_names
  stride <- length(names)
  flat <- c(species$ode_state, species$new_node$ode_state)
  s <- matrix(flat, nrow = stride, dimnames = list(names, NULL))
  data.frame(height = s["height", ],
             area_heartwood = s["area_heartwood", ],
             mass_heartwood = s["mass_heartwood", ],
             log_density = s["log_density", ],
             birth_date = c(vapply(species$nodes,
                                   function(n) n$introduction_time, 0),
                            species$new_node$introduction_time))
}

trapezium_r <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1))) / 2
}

# The census in R. A census is a quadrature of a density, so the grid is the
# coordinate the density is carried in: the birth date, which ascends as the
# nodes are stored and puts the boundary node last, or the height, which ascends
# the other way and puts it first.
census_r <- function(species, pars, eta, include_boundary = TRUE,
                     birth_date = TRUE) {
  st <- species_state_r(species)
  n <- nrow(st)
  order <- if (birth_date) {
    if (include_boundary) seq_len(n) else seq_len(n - 1)
  } else {
    if (include_boundary) c(n, seq(n - 1, 1)) else seq(n - 1, 1)
  }
  st <- st[order, , drop = FALSE]
  psi <- tf24_allometry_r(pars, eta)(st$height, st$area_heartwood,
                                     st$mass_heartwood)
  density <- exp(st$log_density)
  x <- if (birth_date) st$birth_date else st$height
  vapply(psi, function(p) trapezium_r(x, density * p), numeric(1))
}

# A solved stand on the birth-date coordinate, named rather than defaulted: the
# reverse sweep transposes only that one, and it is the grid census_r's reference
# reduction is written on.
solved_stand <- function(lifetime = 20, schedule = NULL) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  # A reverse sweep costs one recording and one sweep per cohort per stage per
  # step, so its cost is set by the NODE COUNT and not by the lifetime. The
  # default schedule gives 81 nodes at a lifetime of 2, where the forward run
  # takes 6 seconds and the sweep takes 1006. A test that needs a realistic size
  # distribution takes the default; one that only needs an entry point to answer
  # passes a short schedule and gets the same code paths for seconds.
  if (!is.null(schedule)) {
    p$node_schedule_times <- schedule
  }
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"),
                                 Control(node_density_in_birth_date = TRUE))
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

test_that("G3: on the birth-date coordinate the weights carry no derivative", {
  # A height is state and a birth date is not, so which one the grid is built
  # from decides whether a cohort's height moves the quadrature as well as the
  # integrand. On this coordinate it does not, and the two differences below
  # are the same number: the weight term the height grid would carry is larger
  # than the whole derivative and of the opposite sign, so a seed built on the
  # wrong grid is not a small error.
  scm <- solved_stand(12)
  species <- scm$patch$species[[1]]
  strategy <- scm$parameters$strategies[[1]]
  st <- species_state_r(species)
  psi_at <- tf24_allometry_r(strategy$pars, strategy$pars$eta)

  k <- floor(nrow(st) / 4)
  eps <- 1e-6
  # Rows run in storage order with the boundary node last, which ascends in
  # birth date and descends in height; negating the height integral is what
  # takes it up its own axis instead.
  leaf_area_at <- function(h, coordinate) {
    hh <- st$height
    hh[k] <- h
    p <- psi_at(hh, st$area_heartwood, st$mass_heartwood)
    y <- exp(st$log_density) * p$leaf_area
    if (coordinate == "birth_date") {
      trapezium_r(st$birth_date, y)
    } else {
      -trapezium_r(hh, y)
    }
  }
  d_dh <- function(coordinate) {
    (leaf_area_at(st$height[k] + eps, coordinate) -
       leaf_area_at(st$height[k] - eps, coordinate)) / (2 * eps)
  }
  over_birth_date <- d_dh("birth_date")
  over_height <- d_dh("height")

  # Non-vacuity: the two grids must disagree, or this proves nothing. Measured
  # here they disagree by a factor of about 16, so a seed built on the wrong one
  # is not a small error.
  expect_true(abs(over_height - over_birth_date) >
                2 * min(abs(over_height), abs(over_birth_date)))

  # The seed the reverse pass is given is the derivative on the coordinate the
  # density is carried in.
  seed <- try(stand_census_state_adjoint(scm), silent = TRUE)
  skip_if(inherits(seed, "try-error"),
          "the census seed is refused at the leaf boundary")
  stride <- length(scm$patch$species[[1]]$new_node$ode_names)
  col <- (k - 1) * stride + 1
  expect_equal(unname(seed["leaf_area", col]), over_birth_date,
               tolerance = 1e-4)
})

test_that("G4: the seed reaches every state a metric reads", {
  # An exact zero in this design is the signature of a missing accumulator and
  # never of true insensitivity, so a state a metric demonstrably reads must have
  # a seed. Two of them do not.
  #
  # `area_stem` sums bark, sapwood and heartwood AREA, and `mass_above_ground`
  # sums the three masses and heartwood MASS, so both read a heartwood state
  # directly and linearly: the seed is the quadrature weight times the density,
  # with no allometry in between, which is what makes this checkable by hand.
  # G1 shows the census VALUE carries them, so the reduction is right and it is
  # the recording of it that is not.
  # A lifetime where the bottom of the distribution is still alive: where a
  # density has underflowed to zero its seed is legitimately zero and the
  # comparison below says nothing.
  scm <- solved_stand(5)
  seed <- stand_census_state_adjoint(scm)
  species <- scm$patch$species[[1]]
  names_i <- species$nodes[[1]]$ode_names
  stride <- length(names_i)

  # Birth-date trapezium weights, boundary node last, as the census integrates.
  b <- c(vapply(species$nodes, function(n) n$introduction_time, 0),
         species$new_node$introduction_time)
  dens <- c(vapply(species$nodes, function(n) exp(n$log_density), 0),
            exp(species$new_node$log_density))
  w <- vapply(seq_along(b), function(i) {
    lo <- if (i > 1) (b[i] - b[i - 1]) / 2 else 0
    hi <- if (i < length(b)) (b[i + 1] - b[i]) / 2 else 0
    lo + hi
  }, 0)
  expected <- (w * dens)[seq_along(species$nodes)]

  live <- expected > 0
  expect_true(any(live))
  for (pair in list(c("area_stem", "area_heartwood"),
                    c("mass_above_ground", "mass_heartwood"))) {
    cols <- (seq_along(species$nodes) - 1) * stride + match(pair[2], names_i)
    expect_equal(unname(seed[pair[1], cols][live]), expected[live],
                 tolerance = 1e-8, info = paste(pair[1], "reads", pair[2]))
  }

  # Non-vacuity: the states that DO have a seed, so a wholesale failure of the
  # recording would not pass this by looking like the defect above.
  h_cols <- (seq_along(species$nodes) - 1) * stride + match("height", names_i)
  expect_true(any(seed["leaf_area", h_cols] != 0))
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
  seed <- try(stand_census_state_adjoint(scm), silent = TRUE)
  skip_if(inherits(seed, "try-error"),
          "the census seed is refused at the leaf boundary")
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

test_that("the Control a gradient is taken at is the entries that move it", {
  # Four move the TRAJECTORY the gradient is taken along. The fifth moves no
  # forward number at all and still decides which rows exist, by refusing a
  # collar response the profit curvature is too small to support -- so two
  # gradients taken at different floors are gradients of different functions for a
  # different reason, and both reasons belong in the same comparison.
  scm <- solved_stand()
  expect_equal(names(gradient_control(scm)),
               c("GSS_tol_abs", "ci_abs_tol", "node_gradient_eps",
                 "schedule_eps", "gradient_curvature_floor"))
  ctrl <- Control()
  expect_equal(unname(gradient_control(scm)),
               c(ctrl$GSS_tol_abs, ctrl$ci_abs_tol, ctrl$node_gradient_eps,
                 ctrl$schedule_eps, ctrl$gradient_curvature_floor))
})

test_that("the trait gradient entry point is reachable", {
  # The symbol exists and stand_gradient reaches it. It cannot return a gradient
  # on a run whose ODE state widens at a node introduction: the reverse sweep
  # carries one lambda of one width and nothing narrows the system to meet an
  # earlier record, so such a run is refused by name rather than swept at a
  # width its records do not have.
  # Two nodes, because this test asks whether the entry point answers and with
  # what names -- not what the numbers are. On the default schedule the same
  # assertions cost a sweep over eighty-one cohorts.
  scm <- solved_stand(5, schedule = list(c(0, 0.63)))
  expect_true(is.function(census_trait_gradient_tf24))
  # A column is named for its species as well as its parameter: a bare name would
  # resolve to species one's column silently on a multi-species stand.
  got <- tryCatch(stand_gradient(scm, traits = "1.lma"), error = identity)
  if (inherits(got, "error")) {
    # Either refusal is by name: the sweep cannot cross an introduction, or the
    # leaf supplies no rows for the output the water channel runs through.
    expect_match(conditionMessage(got),
                 "widens the ODE state|per-layer uptake")
  } else {
    expect_equal(rownames(got$gradient), census_metric_names_tf24())
    expect_equal(colnames(got$gradient), "1.lma")
  }
  # And the bare name refuses, naming the convention rather than only failing.
  expect_error(stand_gradient(scm, traits = "lma"), "species index")
})
