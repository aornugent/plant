# Stand census and its direct (fixed-state) sensitivity term.
#
# The oracle is a finite difference ON A FROZEN STATE: hold the size
# distribution fixed, perturb the strategy parameter, recompute the census,
# difference. That is exact algebra differenced against the same algebra, so it
# should agree at the arithmetic floor rather than at some trajectory-comparison
# tolerance.
#
# Two species throughout. A single-species suite cannot see the column-naming
# defect this checks for.


census_test_state <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"),
                      birth_rate = c(20, 20))
  res <- run_scm(p, ctrl = control(), collect = TRUE)
  list(state = stand_census_state(res), strategies = res$p$strategies)
}

# Trapezium weights on a height-sorted grid, in input order. Independent of the
# C++ implementation; used to check the census value itself.
trapezium_weights_r <- function(h) {
  n <- length(h)
  o <- order(h)
  dh <- diff(h[o])
  ws <- numeric(n)
  ws[-n] <- ws[-n] + dh / 2
  ws[-1] <- ws[-1] + dh / 2
  w <- numeric(n)
  w[o] <- ws
  w
}

perturb <- function(s, par, value) {
  pars <- s$pars
  pars[[par]] <- value
  s$pars <- pars
  s
}

test_that("census value matches an independent quadrature of the allometry", {
  tt <- census_test_state()
  cen <- stand_census(tt$state, tt$strategies)

  # Rebuild the same three totals from the existing allometry expansion.
  by_hand <- c(area_leaf = 0, area_stem = 0, mass_above_ground = 0)
  for (i in seq_along(tt$state)) {
    st <- tt$state[[i]]
    al <- FF16_strategy_expand_allometry(tt$strategies[[i]], st$height,
                                         st$area_heartwood, st$mass_heartwood)
    w <- trapezium_weights_r(st$height) * st$density
    by_hand["area_leaf"] <- by_hand["area_leaf"] + sum(w * al$area_leaf)
    by_hand["area_stem"] <- by_hand["area_stem"] + sum(w * al$area_stem)
    by_hand["mass_above_ground"] <-
      by_hand["mass_above_ground"] + sum(w * al$mass_above_ground)
  }
  expect_equal(as.numeric(cen), as.numeric(by_hand[names(cen)]),
               tolerance = 1e-12)
  expect_true(all(attr(cen, "grid_is_monotone")))
})

test_that("direct term agrees with a frozen-state finite difference", {
  tt <- census_test_state()
  out <- stand_census_direct_term(tt$state, tt$strategies)
  pars <- census_parameters()

  # Central difference at two step sizes, Richardson-extrapolated to remove the
  # O(h^2) truncation term. The binding error here is cancellation, not
  # truncation: above-ground mass is dominated by heartwood, which carries no
  # parameter at all, so d/d(rho) is ~1e-7 of a total of ~14 and a small step
  # loses most of the significant digits of the difference. Killing the h^2
  # term lets the step be large enough for that to stop mattering. The two
  # allometric constants enter through a power law and keep a small step;
  # everything else enters linearly, where a large step costs nothing.
  central <- function(i, par, h) {
    base <- tt$strategies[[i]]$pars[[par]]
    up <- tt$strategies
    dn <- tt$strategies
    up[[i]] <- perturb(up[[i]], par, base + h)
    dn[[i]] <- perturb(dn[[i]], par, base - h)
    as.numeric(stand_census(tt$state, up) - stand_census(tt$state, dn)) / (2 * h)
  }
  step <- function(par) if (par %in% c("a_l1", "a_l2")) 1e-3 else 1e-2

  worst <- setNames(numeric(length(census_metrics())), census_metrics())
  for (i in seq_along(tt$state)) {
    for (par in pars) {
      h <- step(par) * abs(tt$strategies[[i]]$pars[[par]])
      fd <- (4 * central(i, par, h / 2) - central(i, par, h)) / 3
      an <- as.numeric(out$gradient[, paste0(names(tt$state)[[i]], ".", par)])
      scale <- pmax(abs(fd), abs(an))
      worst <- pmax(worst, ifelse(scale > 0, abs(fd - an) / scale, 0))
    }
  }
  # Exact algebra differenced against the same algebra: this is the arithmetic
  # floor, not a modelling tolerance. `area_leaf` has no parameter-free part and
  # so has no cancellation; `mass_above_ground` is where the heartwood constant
  # sets the floor.
  message("census direct term, worst relative FD disagreement: ",
          paste(names(worst), format(worst, digits = 3), collapse = ", "))
  expect_lt(max(worst), 1e-8)
  expect_lt(worst[["area_leaf"]], 1e-10)
})

test_that("a structural zero is marked, and a marked term is not zero", {
  tt <- census_test_state()
  out <- stand_census_direct_term(tt$state, tt$strategies)

  # An absence: reported as an exact zero AND flagged, so it cannot be confused
  # with a term that ran.
  expect_true(all(out$gradient[!out$support] == 0))
  expect_false(all(out$support))

  # Leaf area reads only the two allometric constants.
  la <- out$support["area_leaf", grepl("^species_1\\.", colnames(out$support))]
  expect_equal(names(which(la)), c("species_1.a_l1", "species_1.a_l2"))
  # Stem area reaches the stem-area constants through sapwood and bark.
  sa <- out$support["area_stem", grepl("^species_1\\.", colnames(out$support))]
  expect_equal(sort(names(which(sa))),
               sort(c("species_1.a_l1", "species_1.a_l2",
                      "species_1.theta", "species_1.a_b1")))
  # Above-ground mass reads all seven.
  expect_true(all(out$support["mass_above_ground", ]))

  # In this design an exact zero is otherwise the signature of a missing
  # accumulator, so no supported entry may be exactly zero on a live stand.
  expect_true(all(out$gradient[out$support] != 0))
})

test_that("an unknown parameter is refused by name", {
  tt <- census_test_state()
  expect_error(
    stand_census_direct_term(tt$state, tt$strategies,
                             parameters = c("lma", "wood_density")),
    "wood_density")
  expect_silent(stand_census_direct_term(tt$state, tt$strategies,
                                         parameters = c("lma", "rho")))
})

test_that("columns are named per species", {
  tt <- census_test_state()
  out <- stand_census_direct_term(tt$state, tt$strategies)
  cn <- colnames(out$gradient)

  expect_equal(length(cn), length(tt$state) * length(census_parameters()))
  expect_equal(anyDuplicated(cn), 0L)
  expect_true(all(c("species_1.lma", "species_2.lma") %in% cn))

  # The two species differ in lma, so their columns must differ. Without the
  # species prefix both columns would be called "lma" and `[, "lma"]` would
  # resolve to species one for both.
  expect_false(isTRUE(all.equal(out$gradient[, "species_1.lma"],
                                out$gradient[, "species_2.lma"])))
})

test_that("the census quadrature is guarded against a crossed grid", {
  tt <- census_test_state()
  state <- tt$state
  strategies <- tt$strategies

  st <- state[[1]]
  n <- nrow(st)

  # The trapezium taken in node order, as an unguarded quadrature would take it.
  # On a crossed grid neighbouring trapezia have opposite-signed widths.
  unguarded_weights <- function(h) {
    m <- length(h)
    dh <- diff(h)
    w <- numeric(m)
    w[-m] <- w[-m] - dh / 2
    w[-1] <- w[-1] - dh / 2
    w
  }
  leaf_area_two_ways <- function(sd) {
    al <- FF16_strategy_expand_allometry(strategies[[1]], sd$height,
                                         sd$area_heartwood, sd$mass_heartwood)
    c(unguarded = sum(unguarded_weights(sd$height) * sd$density * al$area_leaf),
      guarded = sum(trapezium_weights_r(sd$height) * sd$density * al$area_leaf))
  }

  # A younger node overtakes an older one: swap two adjacent heights, leaving
  # each node's density and heartwood attached to it. This is the shape
  # reserve-gated growth produces (#517, #571). Do it among the shortest,
  # most abundant cohorts, which is where the census actually has its mass and
  # where the reserve gate bites.
  crossed <- st
  crossed$height[c(n - 2, n - 1)] <- st$height[c(n - 1, n - 2)]
  state[[1]] <- crossed

  guarded <- stand_census(state, strategies)
  expect_false(all(attr(guarded, "grid_is_monotone")))

  two <- leaf_area_two_ways(crossed)
  err <- abs(two[["unguarded"]] - two[["guarded"]]) / abs(two[["guarded"]])
  message("crossed-grid leaf-area census error without the guard: ",
          format(100 * err, digits = 3), "%")
  expect_gt(err, 0.05)

  # One node further down the grid the cancellation is larger than the census:
  # the unguarded quadrature reports a negative total leaf area.
  worse <- st
  worse$height[c(n - 1, n)] <- st$height[c(n, n - 1)]
  two_worse <- leaf_area_two_ways(worse)
  message("one node lower, unguarded leaf area = ",
          format(two_worse[["unguarded"]], digits = 3), " vs ",
          format(two_worse[["guarded"]], digits = 3))
  expect_lt(two_worse[["unguarded"]], 0)
  expect_gt(two_worse[["guarded"]], 0)

  # The guarded census sees only the state, not the node order.
  sorted <- state
  sorted[[1]] <- crossed[order(crossed$height), , drop = FALSE]
  expect_equal(as.numeric(stand_census(sorted, strategies)),
               as.numeric(guarded), tolerance = 1e-14)
})

test_that("the state really is frozen under a parameter perturbation", {
  tt <- census_test_state()
  before <- tt$state
  strategies <- tt$strategies

  # Weights and densities computed once, from the state alone.
  w <- lapply(before, function(st) trapezium_weights_r(st$height) * st$density)

  bumped <- strategies
  bumped[[1]] <- perturb(bumped[[1]], "lma", strategies[[1]]$pars$lma * 1.1)

  base <- stand_census(before, strategies)
  pert <- stand_census(before, bumped)

  # Nothing the perturbation did touched the grid.
  expect_identical(before, tt$state)
  for (i in seq_along(before)) {
    expect_identical(w[[i]],
                     trapezium_weights_r(before[[i]]$height) * before[[i]]$density)
  }

  # And the whole change in the census is the change in the per-individual
  # algebra, re-weighted by those same frozen weights.
  delta <- 0
  for (i in seq_along(before)) {
    st <- before[[i]]
    a <- FF16_strategy_expand_allometry(strategies[[i]], st$height,
                                        st$area_heartwood, st$mass_heartwood)
    b <- FF16_strategy_expand_allometry(bumped[[i]], st$height,
                                        st$area_heartwood, st$mass_heartwood)
    delta <- delta + sum(w[[i]] * (b$mass_above_ground - a$mass_above_ground))
  }
  # 1e-9 rather than machine epsilon because both sides are a difference of two
  # ~14 kg m-2 totals: the change is 3.3e-5 of the census, so the difference
  # keeps ~11 digits, not 16.
  expect_equal(as.numeric(pert[["mass_above_ground"]]) -
                 as.numeric(base[["mass_above_ground"]]),
               delta, tolerance = 1e-9)
})
