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

    cmp    <- local_error_integration(sp$heights, cmp_competition_effect, 1.0)
    cmp_pi <- local_error_integration(sp$heights, cmp_competition_effect, pi)

    expect_identical(sp$compute_competition_effect_by_nodes_error(1.0), cmp)
    expect_identical(sp$compute_competition_effect_by_nodes_error(pi), cmp_pi)

    ode_size <- Node(x, e)(strategy_types[[x]]())$ode_size
    ode_state <- sp$ode_state
    nodes <- sp$nodes
    expect_equal(length(ode_state), ode_size * sp$size)
    expect_identical(ode_state, unlist(lapply(nodes, function(p) p$ode_state)))
  })

  ## The transport term of the size-density equation, on the cohort grid.
  node_ode <- function() Node(x, e)(strategy_types[[x]]())
  three_nodes <- function() {
    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    for (i in seq_len(3)) sp$introduce_new_node()
    sp$heights <- sp$height_max * 40 * c(1, .8, .6)
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    list(sp = sp, env = env)
  }
  growth_rates <- function(sp) {
    vapply(seq_len(sp$size),
           function(i) sp$node_at(i)$individual$rate("height"), numeric(1))
  }

  test_that("the transport term differences growth across neighbouring nodes", {
    f <- three_nodes()
    sp <- f$sp
    g <- growth_rates(sp)
    h <- sp$heights

    expect_equal(sp$growth_rate_gradient(1), (g[1] - g[2]) / (h[1] - h[2]))
    expect_equal(sp$growth_rate_gradient(2), (g[2] - g[3]) / (h[2] - h[3]))

    ## The lowest node differences against the inflow boundary node.
    expect_equal(sp$growth_rate_gradient(3),
                 (g[3] - sp$new_node$individual$rate("height")) /
                   (h[3] - sp$new_node$height))

    ## And it is what the ODE rate carries, less mortality.
    n <- node_ode()$ode_size
    slot <- which(node_ode()$ode_names == "log_density")
    ldd <- sp$ode_rates[(seq_len(sp$size) - 1) * n + slot]
    mort <- vapply(seq_len(sp$size),
                   function(i) sp$node_at(i)$individual$rate("mortality"),
                   numeric(1))
    expect_equal(ldd,
                 -vapply(seq_len(sp$size), sp$growth_rate_gradient, numeric(1)) -
                   mort)

    expect_error(sp$growth_rate_gradient(4), "out of bounds")
  })

  test_that("a node with no interval below it takes the value above", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)

    ## Every node introduced and none moved: all three sit at height_0, as does
    ## the boundary node, so no node has an interval below it.
    for (i in seq_len(3)) sp$introduce_new_node()
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    expect_identical(sp$growth_rate_gradient(1), 0.0)
    expect_identical(sp$growth_rate_gradient(2), 0.0)
    expect_identical(sp$growth_rate_gradient(3), 0.0)

    ## One coincident pair among nodes that are otherwise apart.
    f <- three_nodes()
    sp <- f$sp
    sp$heights <- sp$height_max * c(1, .8, .8)
    sp$compute_rates(f$env, pr_patch_survival = 1, birth_rate = 1)
    expect_identical(sp$growth_rate_gradient(2), sp$growth_rate_gradient(1))

    ## Heights need not stay descending: the ODE state is written without an
    ## order check, and then list order is not height order.
    f <- three_nodes()
    sp <- f$sp
    n <- node_ode()$ode_size
    slot <- which(node_ode()$ode_names == "height")
    y <- sp$ode_state
    y[2 * n + slot] <- y[1 * n + slot] + 1
    sp$ode_state <- y
    sp$compute_rates(f$env, pr_patch_survival = 1, birth_rate = 1)
    expect_lt(sp$heights[2] - sp$heights[3], 0)
    expect_identical(sp$growth_rate_gradient(2), sp$growth_rate_gradient(1))
  })

  test_that("the transport term is the compression rate of the interval below", {
    f <- three_nodes()
    sp <- f$sp
    n <- node_ode()$ode_size
    slot <- which(node_ode()$ode_names == "log_density")

    rates <- sp$ode_rates
    ldd <- rates[(seq_len(sp$size) - 1) * n + slot]
    mort <- vapply(seq_len(sp$size),
                   function(i) sp$node_at(i)$individual$rate("mortality"),
                   numeric(1))
    before <- sp$heights

    dt <- 1e-6
    sp$ode_state <- sp$ode_state + dt * rates
    after <- sp$heights

    ## Interior intervals only: both of their endpoints are transported heights.
    ## The boundary node stays at height_0 and is not one.
    interior <- seq_len(sp$size - 1)
    dlog_dh <- (log(after[interior] - after[interior + 1]) -
                log(before[interior] - before[interior + 1])) / dt
    expect_equal(ldd[interior] + mort[interior], -dlog_dh, tolerance = 1e-5)
  })

  test_that("set_birth_state restores per-node birth bookkeeping", {
    env <- Environment(x)
    sp <- Species(x, e)(strategy_types[[x]]())
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    for (i in seq_len(3)) sp$introduce_new_node()
    sp$heights <- sp$height_max * 40 * c(1, .75, .6)

    times   <- c(0.25, 0.50, 0.75)
    density <- c(1.5, 2.5, 3.5)
    pr_surv <- c(0.9, 0.8, 0.7)
    sp$set_birth_state(times, density, pr_surv)

    ## Each argument reaches its own destination, in node order: the three are
    ## same-typed and unrelated in meaning, so a swap would only show here.
    expect_identical(sp$node_times, times)
    expect_identical(sp$patch_densities, density)
    expect_identical(sp$pr_patch_survival_at_birth, pr_surv)

    ## None of the three is part of the ODE state, so setting them moves none of it.
    state <- sp$ode_state
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    rates <- sp$ode_rates
    sp$set_birth_state(times * 2, density * 2, pr_surv)
    expect_identical(sp$ode_state, state)
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    expect_identical(sp$ode_rates, rates)

    ## pr_patch_survival_at_birth divides the fecundity rate, which is the one
    ## rate it reaches: halving it doubles that rate and leaves the rest alone.
    sp$set_birth_state(times, density, pr_surv / 2)
    sp$compute_rates(env, pr_patch_survival = 1, birth_rate = 1)
    halved <- sp$ode_rates
    n <- Node(x, e)(strategy_types[[x]]())$ode_size
    fecundity <- (seq_len(sp$size) - 1) * n + (n - 1)
    expect_true(all(rates[fecundity] > 0))
    expect_identical(halved[fecundity], rates[fecundity] * 2)
    expect_identical(halved[-fecundity], rates[-fecundity])

    ## Each length must match the node count.
    expect_error(sp$set_birth_state(times[-1], density, pr_surv), "Incorrect length input")
    expect_error(sp$set_birth_state(times, density[-1], pr_surv), "Incorrect length input")
    expect_error(sp$set_birth_state(times, density, pr_surv[-1]), "Incorrect length input")
  })
}
