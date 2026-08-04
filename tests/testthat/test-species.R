strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  e <- environment_types[[x]]


  test_that("Basics", {
    env <- Environment(x)
    s <- strategy_types[[x]]()
    sp <- Species(x, e)(s)
    new_node <- Node(x, e)(s)
    plant <- Individual(x, e)(s)
    h0 <- new_node$height

    expect_equal(sp$size, 0)
    expect_identical(sp$height_max, h0)
    expect_identical(sp$nodes, list())
    expect_identical(sp$height, NULL)
    expect_identical(sp$log_densities, numeric(0))
    expect_identical(sp$compute_competition_effect_by_nodes, numeric(0))
    expect_identical(sp$compute_competition_effect_by_nodes_error(1.0), numeric(0))
    expect_equal(sp$ode_size, 0)
    expect_identical(sp$ode_state, numeric(0))
    expect_identical(sp$ode_rates, numeric(0))

    ## Causes initial conditions to be estimated:
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    new_node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)

    ## Internal and test new_node report same values:
    expect_identical(sp$new_node$rates, new_node$rates)
    expect_identical(sp$new_node$ode_state, new_node$ode_state)

    sp$introduce_new_node()
    expect_equal(sp$size, 1)

    nodes <- sp$nodes
    expect_inherits(nodes, "list")
    expect_equal(length(nodes), 1)
    expect_identical(nodes[[1]]$rates, new_node$rates)
    expect_equal(sp$heights, new_node$height)
    expect_equal(sp$log_densities, new_node$log_density)
    expect_equal(sp$compute_competition_effect_by_nodes, new_node$compute_competition(0.0))
    ## NOTE: Didn't check ode values

    ## Internal and test new_node report same values:
    expect_identical(sp$new_node$rates, new_node$rates)

    expect_inherits(sp$node_at(1), sprintf("Node<%s,%s>",x,e))
    expect_identical(sp$node_at(1)$rates, nodes[[1]]$rates)

    ## Not sure about this -- do we need more immediate access?
    expect_identical(sp$new_node$individual$establishment_probability(env), plant$establishment_probability(env))

    expect_equal(sp$compute_competition(0), 0)

    sp$heights <- 1

    h <- 0
    xx <- c(sp$new_node$height, sp$heights)
    y <- c(sp$new_node$compute_competition(h),
           sp$node_at(1)$compute_competition(h))

    expect_identical(sp$compute_competition(h), trapezium(xx, y))

    ## Better tests: I want cases where:
    ## 1. empty: throws error
    sp$clear()

    ## Re-set up the initial conditions
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)

    expect_error(sp$node_at(1), "Index 1 out of bounds")
    expect_error(sp$node_at(0), "Invalid value for index")
  })

  ## 1: empty species (no nodes) has no leaf area above any height:
  test_that("Empty species has no leaf area", {
    sp <- Species(x, e)(strategy_types[[x]]())
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  ## 2: Node up against boundary has no leaf area:
  test_that("species with only boundary node no leaf area", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$introduce_new_node()
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    expect_equal(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(10), 0)
    expect_equal(sp$compute_competition(Inf), 0)
  })

  cmp_compute_competition <- function(h, sp) {
    x <- c(sp$heights, sp$new_node$height)
    y <- c(sapply(sp$nodes, function(p) p$compute_competition(h)),
           sp$new_node$compute_competition(h))
    trapezium(rev(x), rev(y))
  }

  ## The local trapezium error over the grid compute_competition() actually
  ## integrates: the node abscissae *closed by the boundary node's*, which is
  ## where its last trapezium ends (see cmp_compute_competition above). The
  ## metric drives schedule refinement, so it has to be taken over that grid and
  ## not over the node list alone.
  cmp_competition_error <- function(sp, scal) {
    x <- c(sp$heights, sp$new_node$height)
    y <- c(sapply(sp$nodes, function(p) p$compute_competition(0.0)),
           sp$new_node$compute_competition(0.0))
    local_error_integration(x, y, scal)
  }

  ## 3: Single node; one round of trapezium:
  test_that("Leaf area sensible with one node", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    sp$introduce_new_node()
    h_top <- sp$height_max * 4
    sp$heights <- h_top

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))

    expect_identical(sp$compute_competition(h_top), 0.0)

    ## Part way up (and above bottom offspring boundary condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))

    ode_size <- Node(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    p <- sp$node_at(1)
    expect_equal(sp$ode_size, ode_size)
    expect_equal(length(ode_state), ode_size)
    expect_identical(ode_state, p$ode_state)
  })

  test_that("Leaf area sensible with two nodes", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    sp$introduce_new_node()
    h_top <- sp$height_max * 4
    sp$introduce_new_node()
    sp$heights <- h_top * c(1, .6)

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom node, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the offspring)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    ode_size <- Node(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    nodes <- sp$nodes
    expect_equal(sp$ode_size, ode_size * sp$size)
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(nodes, function(p) p$ode_state)))
  })

  test_that("Leaf area sensible with three nodes", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    sp$introduce_new_node()
    h_top <- sp$height_max * 4
    sp$introduce_new_node()
    sp$introduce_new_node()
    sp$heights <- h_top * c(1, .75, .6)

    ## At base and top
    expect_gt(sp$compute_competition(0), 0)
    expect_equal(sp$compute_competition(0), cmp_compute_competition(0, sp))
    expect_equal(sp$compute_competition(h_top), 0)
    ## Part way up (below bottom node, above boundarty condition)
    expect_equal(sp$compute_competition(h_top * .5), cmp_compute_competition(h_top * .5, sp))
    ## Within the top pair (excluding the offspring)
    expect_equal(sp$compute_competition(h_top * .8), cmp_compute_competition(h_top * .8, sp))

    cmp_competition_effect <- sapply(seq_len(sp$size),
                            function(i) sp$node_at(i)$compute_competition(0.0))
    expect_identical(sp$compute_competition_effect_by_nodes, cmp_competition_effect)

    cmp    <- cmp_competition_error(sp, 1.0)
    cmp_pi <- cmp_competition_error(sp, pi)

    expect_identical(sp$compute_competition_effect_by_nodes_error(1.0), cmp)
    expect_identical(sp$compute_competition_effect_by_nodes_error(pi), cmp_pi)

    ## One entry per point of that grid, the boundary node included, so the
    ## youngest node -- the only one whose error depends on the closing
    ## segment -- gets a real error rather than the NA of a grid end.
    expect_equal(length(cmp), sp$size + 1L)
    expect_true(is.finite(cmp[[sp$size]]))

    ode_size <- Node(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    nodes <- sp$nodes
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(nodes, function(p) p$ode_state)))
  })
}

## Species::compute_competition() integrates a trapezium over the node list and
## then closes it with a segment running from the youngest node to the boundary
## node, because the birth boundary -- not the youngest cohort -- is where the
## size distribution starts. That segment is not a rounding detail: density is
## peaked at the boundary, and it carries a few percent of the competition
## integral on a refined schedule and far more early in a run.
##
## The refinement error metric is a local trapezium error over that same
## quadrature, so it has to be measured over the same grid. Over the node list
## alone the boundary node's abscissa is missing, the youngest node's error comes
## back NA, and no error at all is reported for the closing segment -- the
## adaptive refiner cannot see the quadrature it is actually performing. Checked
## in both coordinates, since Species::quadrature_abscissa switches between
## introduction_time() and -height(), and on real patch states, which is where
## the boundary node's abscissa has to be current (Patch::compute_environment
## restamps it) rather than left over from an earlier derivs call.
test_that("the refinement error metric spans the closing segment", {
  x <- "FF16"; e <- "FF16_Env"
  p0 <- scm_base_parameters(x)
  p0$max_patch_lifetime <- 20
  p <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)

  for (birth_date in c(FALSE, TRUE)) {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- birth_date
    scm <- SCM(x, e)(p, Environment(x), ctrl)
    scm$collect <- TRUE
    scm$run()

    shares <- numeric(0)
    errors <- numeric(0)
    for (patch in scm$history) {
      s <- patch$species[[1]]
      n <- s$size
      f <- s$compute_competition_effect_by_nodes
      ## Two points make no local error, and on the height path a youngest node
      ## whose density has collapsed to exactly zero drops the closing segment
      ## from the integral too -- then the node list really is the whole grid.
      if (n < 3L || !(birth_date || f[[n]] > 0)) next

      tot <- s$compute_competition(0) / patch$get_area
      abscissa <- if (birth_date) {
        c(s$node_times, s$new_node$introduction_time)
      } else {
        c(s$heights, s$new_node$height)
      }
      effect <- c(f, s$new_node$compute_competition(0))
      err <- s$compute_competition_effect_by_nodes_error(tot)

      ## One entry per point of the quadrature grid, and the same numbers a
      ## trapezium local error over that grid gives.
      expect_equal(length(err), n + 1L)
      expect_equal(err, local_error_integration(abscissa, effect, tot))
      ## The youngest node's entry is the one the closing segment enters; it is
      ## NA when the boundary node is left out of the grid.
      expect_true(is.finite(err[[n]]))

      shares <- c(shares,
                  abs(diff(tail(abscissa, 2))) * sum(tail(effect, 2)) / (2 * tot))
      errors <- c(errors, err[[n]])
    }

    ## And it earns its place in the metric. The segment carries a large share of
    ## the integral over much of the run (a third of it at the youngest states,
    ## ~10% at the median), and among those states the error newly reported for it
    ## is orders of magnitude above round-off -- measured 1e-3 (height) and 2e-3
    ## (birth date) against interior errors of 9e-3 and 1.4e-2, so the same order
    ## as the errors the metric already reported, not a formality. A max over the
    ## large-share states rather than the single largest-share one, because that
    ## one is the youngest state, where the first cohorts sit almost on top of the
    ## boundary node: widest share, but a near-linear integrand and a local error
    ## of ~1e-10.
    expect_gt(max(shares), 0.05)
    expect_gt(max(errors[shares > 0.05]), 1e-4)
  }
})
