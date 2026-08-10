# Rung 5: introductions, where the state changes dimension.
#
# The minimum fixture is two species and three introductions in the order
# species 1, species 2, species 1 -- four segments, with the node stride
# exercised in both directions. Introducing into species 1 while species 2 exists
# is the arrangement a narrowing implemented as a truncation of the tail fails
# on.
#
# Each boundary channel is probed by a parameter reaching the census through that
# channel and no other, so a dropped channel gives an exact zero rather than a
# small error. That is the cleanest test shape in the ladder, and the boundary is
# where it is available.

test_that("the fixture actually widens, in both stride directions", {
  # Non-vacuity before anything else: a probe of the introduction machinery on a
  # stand that introduces nothing proves nothing.
  stand <- ladder_stand_introductions()
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  changes <- which(diff(widths) != 0)
  message(sprintf("\n  %d recorded steps, %d widenings, widths %s",
                  length(widths), length(changes),
                  paste(unique(widths), collapse = " -> ")))
  expect_gte(length(changes), 3L)

  # Species 1 is introduced into while species 2 already exists, which is what
  # makes a tail truncation wrong rather than merely narrow.
  nodes <- ladder_nodes(stand)
  first <- nodes[nodes$species == 1 & !nodes$boundary, ]
  second <- nodes[nodes$species == 2 & !nodes$boundary, ]
  expect_gt(max(first$birth_date), min(second$birth_date))
})

test_that("the initial reserve channel carries a row", {
  # A germination event provisions the seedling: the initial reserve is a trait
  # times the storage capacity, so how well a seed is stocked reads the traits
  # directly. It enters nowhere else, so a dropped channel is an exact zero.
  stand <- ladder_stand_introductions()
  result <- ladder_gradient_or_skip(stand)
  at <- which(census_trait_names_tf24(stand) == "a_st3")
  expect_gt(length(at), 0L)
  peak <- max(abs(result$gradient[, at]))
  message(sprintf("\n  a_st3 row peak: %.6e", peak))
  expect_gt(peak, 0)
})

test_that("the newcomer depends on the state it was introduced into", {
  # The introduction is evaluated at the pre-introduction state: the field is
  # rebuilt and the boundary node placed in it, so the newcomer's own quantities
  # depend on the state before the introduction as well as on the traits. An
  # implementation written as though the newcomer has no predecessor drops that
  # half entirely.
  #
  # Perturbing the soil at the introduction time is what exposes it, because a
  # boundary node taking its potentials from a cache has no visible dependence on
  # moisture at all -- and establishment is the most moisture-sensitive event in
  # the life cycle.
  base <- ladder_stand_introductions()
  ladder_gradient_or_skip(base)

  perturbed <- local({
    p <- ladder_parameters(c("fast", "slow"))
    p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
    env <- Environment("TF24")
    env$set_soil_water_state(rep(0.26, env$get_soil_number_of_depths()))
    run_scm(p, env, ladder_control(), collect = FALSE)
  })
  moved <- stand_gradient(perturbed)

  ladder_expect_moves(stand_gradient(base)$gradient, moved$gradient,
                      "the newcomer's dependence on pre-introduction state")
})

test_that("the initial condition's own adjoint is not discarded", {
  # The last segment's state adjoint is the initial-condition term. A sweep that
  # narrows through every introduction and then discards the final narrowed
  # vector has computed it and thrown it away, and the reserve channel has no
  # other path.
  #
  # The probe is a parameter reaching only the founding cohorts: with
  # introductions after time zero, a row that exists for a later cohort and not
  # for a founding one is the signature.
  stand <- ladder_stand_introductions()
  with_founders <- ladder_gradient_or_skip(stand)

  later_only <- local({
    p <- ladder_parameters(c("fast", "slow"))
    p$node_schedule_times <- list(c(0.29, 0.94), c(0.57))
    ladder_run(p)
  })
  moved <- stand_gradient(later_only)
  ladder_expect_moves(with_founders$gradient, moved$gradient,
                      "the founding cohorts' contribution")
})

test_that("seed height is on a declared list rather than silently zero", {
  # The seed height solves an implicit condition on the strategy and the
  # strategy's preparation is evaluated in double, so every parameter reaching
  # birth size carries no row through that channel. This is imposed, not derived,
  # and it is the one term no available instrument can referee: a tangent imposes
  # the same equation and a re-run finite difference cannot referee at production
  # because a relative step of two parts in ten million moves a mature stand
  # between alive and identically zero.
  #
  # So it is not a numerical probe. What it must be is declared, because a silent
  # zero here is the failure. Eight parameters reach birth size, and each keeps
  # whatever rows its other routes give it, so this is a declared bias rather
  # than a column of zeros.
  declared <- ladder_birth_size_channel_zero()
  expect_length(declared, 8L)

  stand <- ladder_stand_introductions()
  columns <- census_trait_names_tf24(stand)
  expect_length(setdiff(declared, columns), 0L)
})

test_that("a marginal recruit's gradient exists, is finite, and tends to zero", {
  # Both the cumulative mortality and the log density diverge logarithmically as
  # establishment goes to zero, and the sensitivity of the log density diverges
  # with them. But the log density reaches everything downstream only through the
  # stem number, whose own sensitivity tends to zero, and the census carries stem
  # number linearly. So the finite answer is reached by seeding the quantity that
  # has one.
  #
  # The fixture assertion is derived rather than chosen: the establishment
  # probability's derivative peaks at a known point, and a recruit anywhere else
  # is in a flat region where the check proves nothing.
  skip_if_not(
    exists("ladder_establishment_probe_tf24", mode = "function"),
    paste("nothing reports a boundary node's establishment probability or the",
          "production driving it, so a recruit cannot be placed in the stiff",
          "band and this check cannot be made non-vacuous"))
})
