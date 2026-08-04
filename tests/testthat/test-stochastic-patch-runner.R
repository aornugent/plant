
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

## ODE states contributed by the environment; only TF24 carries any (five
## soil-water layers then four cumulative fluxes).
n_environment_ode_states <- c(FF16 = 0L, TF24 = 9L, K93 = 0L)

test_that("empty", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()))
    
    env <- Environment(x)
    ctrl <- Control()

    obj <- StochasticPatchRunner(x, e)(p, env, ctrl)
    expect_identical(obj$time, 0.0)

    sched <- obj$node_schedule
    expect_equal(sched$size, 0)
    expect_equal(sched$max_time, p$max_patch_lifetime)

    ## Now, create a new set of times:
    sched2 <- stochastic_schedule(p)
    expect_gt(sched2$size, 0)

    ## Does this need to happen twice?
    obj$node_schedule <- sched2
    expect_equal(obj$node_schedule$size, sched2$size)

    ## Importantly, this moves time forward to where the first
    ## introduction will be!
    expect_identical(obj$time, sched2$next_event$time_introduction)

    ## We're empty though.... The ODE system still carries the environment's own
    ## states (TF24's soil water and cumulative fluxes; none for FF16 or K93).
    expect_equal(obj$patch$species[[1]]$size, 0)
    expect_equal(obj$patch$node_ode_size, 0)
    expect_length(obj$patch$ode_state, n_environment_ode_states[[x]])

    res <- obj$run_next()
    expect_equal(res, 1L)
    expect_identical(obj$time, sched2$all_times[[1]][[2]])

    ode_size <- Individual(x, e)(strategy_types[[x]]())$ode_size
    expect_equal(obj$patch$node_ode_size, ode_size)
    expect_equal(length(obj$patch$ode_state),
                 ode_size + n_environment_ode_states[[x]])
    expect_equal(obj$patch$size, 1)

    expect_false(obj$complete)
  }
})

test_that("collect returns a well-formed, non-empty trajectory (#498)", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    set.seed(1)
    p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                          patch_area=50)
    res <- run_stochastic_collect(p, Environment(x), Control())

    ## Regression guard for #498: the collector used to read a removed `state`
    ## accessor and silently returned empty output, which `expect_silent` could
    ## not catch. Assert the trajectory is actually populated.
    expect_setequal(names(res), c("time", "species", "light_env", "p"))
    expect_gt(length(res$time), 1)
    expect_false(is.unsorted(res$time))          # patch age is non-decreasing
    expect_length(res$species, 1)

    sp <- res$species[[1]]
    expect_equal(length(dim(sp)), 3)             # [variable, time, plant]
    expect_equal(dim(sp)[2], length(res$time))
    expect_true("height" %in% dimnames(sp)[[1]])

    ia <- attr(res$species, "is_alive")[[1]]
    expect_equal(dim(ia), c(length(res$time), dim(sp)[3]))  # [time, plant]

    ## At least one individual was introduced and some survive to the end.
    expect_gt(dim(sp)[3], 0)
    expect_gt(sum(ia[nrow(ia), ], na.rm = TRUE), 0)

    ## Heights are finite wherever an individual is alive. Padded (not-yet-born)
    ## cells are NA in `is_alive`, so select with which() to drop them.
    heights <- sp["height", , ]
    expect_true(all(is.finite(heights[which(ia == 1)])))
  }
})

test_that("collect output is reproducible and matches a seeded baseline (#482)", {
  ## With set.seed(1) and patch_area = 50 the arrival schedule (R RNG) and the
  ## deaths (R::unif_rand, in C++) are fully reproducible, so the number of
  ## individuals introduced and the number alive at the final step are fixed.
  ## These golden values guard against trajectory-changing regressions in the
  ## stochastic tower; update them deliberately if the model/RNG use changes.
  ##
  ## n_total is fixed by the R-side arrival draw and so is insensitive to the
  ## model. TF24's n_alive_final moved from 28 to 31 when the patch started
  ## integrating the environment's own ODE state: the number of mortality draws
  ## is unchanged, but their probabilities are not, because soil water now
  ## evolves instead of sitting at its initial value. FF16 and K93 carry no
  ## environment state, so their counts must not move.
  baseline <- list(
    FF16 = list(n_total = 103L, n_alive_final = 25L),
    TF24 = list(n_total = 103L, n_alive_final = 31L),
    K93  = list(n_total = 105L, n_alive_final = 42L)
  )
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    run_once <- function() {
      set.seed(1)
      p <- Parameters(x, e)(strategies=list(strategy_types[[x]]()),
                            patch_area=50)
      run_stochastic_collect(p, Environment(x), Control())
    }
    res <- run_once()
    sp <- res$species[[1]]
    ia <- attr(res$species, "is_alive")[[1]]
    expect_equal(dim(sp)[3], baseline[[x]]$n_total)
    expect_equal(sum(ia[nrow(ia), ], na.rm = TRUE), baseline[[x]]$n_alive_final)

    ## Determinism: same seed -> identical trajectory.
    res2 <- run_once()
    expect_equal(res$time, res2$time)
    expect_equal(res$species, res2$species)
  }
})

test_that("TF24's environment states are integrated over a stochastic run", {
  ## The stochastic patch used to run the species alone, so TF24's soil water sat
  ## at its initial 0.214 for a whole run while the SCM integrated it up to the
  ## drainage equilibrium and then drew it down as leaf area grew.
  p <- Parameters("TF24", "TF24_Env")(strategies=list(TF24_Strategy()),
                                      patch_area=50)
  p$max_patch_lifetime <- 8
  env <- Environment("TF24")
  expect_equal(env$get_soil_water_state(), rep(0.214, 5))

  set.seed(1)
  obj <- StochasticPatchRunner("TF24", "TF24_Env")(p, env, Control())
  obj$set_node_schedule_times(list(seq(0.2, 7.8, by = 0.2)))
  while (!obj$complete) obj$run_next()
  expect_gt(obj$patch$species[[1]]$size, 0)

  theta <- obj$patch$environment$get_soil_water_state()
  expect_true(all(theta > 0.214))
  expect_true(all(theta < obj$patch$environment$soil_moist_sat))

  ## Cumulative rainfall, infiltration, deep drainage and root uptake. A non-zero
  ## uptake is the individuals' consumption reaching the water balance, which the
  ## per-node and per-species consumption_rate forwarders carry; it cannot exceed
  ## what infiltrated.
  flux <- obj$patch$environment$get_soil_water_state_cumulative_flux()
  expect_length(flux, 4)
  expect_true(all(flux > 0))
  expect_lt(flux[[4]], flux[[2]])
})
