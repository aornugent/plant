strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  s <- strategy_types[[x]]()
  e <- environment_types[[x]]

  test_that("setup, growth rates", {

    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    expect_inherits(node, sprintf("Node<%s,%s>", x, e))
    expect_inherits(node$individual, sprintf("Individual<%s,%s>", x, e))

    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)

    ## The node's own rates are the individual's. The transport term of the
    ## size-density equation is not among them: it needs the node's neighbours,
    ## and lives on Species (see test-species.R).
    plant$compute_rates(env)
    node$compute_rates(env, pr_patch_survival = 1)
    expect_identical(node$individual$rate("height"), plant$rate("height"))
    expect_identical(node$individual$rate("mortality"), plant$rate("mortality"))
  })

  ## TODO(#482): Not done yet:
  ##   * Check that the initial conditions are actually correct
  ##   * Check that the rates computed are actually correct
  test_that("ODE interface", {
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)

    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)
    ## Seed the reference plant's strategy-specific initial states the same way
    ## the node does (via compute_initial_conditions), so the two agree for
    ## strategies that carry seeded states -- e.g. TF24's NSC storage pool (#517).
    plant$set_initial_states(env)
    plant$compute_rates(env)

    nms <- c(plant$ode_names,
             "offspring_produced_survival_weighted", "log_density")
    expect_equal(node$ode_size, length(nms))
    expect_equal(node$ode_names, nms)

    ## Mortality is different because that's what Nodes track
    for( v in setdiff(plant$ode_names, "mortality")) {
      expect_equal(node$individual$state(v), plant$state(v))
    }

    ## Set up plant too:
    pr_estab <- plant$establishment_probability(env)

    y <- plant$ode_state
    g <- plant$rate("height")

    ## Ode *values*:
    cmp <- c(plant$internals$states,
             0, # offspring_produced_survival_weighted
             log(pr_estab / g) # log density
             )
    cmp[which(plant$ode_names == 'mortality')] <- -log(pr_estab)
    expect_equal(node$ode_state, cmp)

    expect_identical(node$fecundity, 0.0);

    ## Ode *rates*:
    cmp <- c(plant$internals$rates,
             ## This is different to the approach in tree1?
             plant$rate("fecundity") * exp(-plant$state("mortality")),
             ## A node held outside a Species has no interval below it, so it
             ## has no transport term.
             0)


    expect_equal(node$ode_rates, cmp)
  })

  test_that("leaf area calculations", {
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    env <- Environment(x)
    env$set_fixed_environment(1, 100)

    h <- node$height

    expect_equal(node$log_density, -Inf) # zero
    expect_equal(exp(node$log_density), 0.0) # zero
    expect_equal(node$compute_competition(0.0), 0) # zero density
    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)

    expect_equal(node$ode_state[[node$ode_size]], node$log_density)
    density <- exp(node$log_density)
    expect_equal(node$compute_competition(0.0), plant$compute_competition(0.0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)

    h <- 8.0
    plant$set_state("height", h)
    v <- node$ode_state
    v[[1]] <- h
    node$ode_state <- v
    expect_identical(plant$state("height"), h)
    expect_identical(node$height, h)

    expect_equal(node$compute_competition(0.0), plant$compute_competition(0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)
  })
}
