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

test_that("a reloaded state carries the boundary node the run carries", {
  # A stage evaluates the inflow condition twice, in two different fields, and
  # they are the same function at different arguments: the first is taken with
  # every species' boundary interval left off, the field is then rebuilt including
  # it, and the second is taken in that rebuilt field. The reductions were built
  # on the first; the water aggregation, an introduced node and a census all read
  # the second.
  #
  # Every rebuild the sweep does -- the census seed, the census's direct term, the
  # replay of an introduction, and the introduction's own transpose -- starts from
  # a recorded state, and a state load alone stops at the first evaluation. A sweep
  # that stops there linearises a boundary node the trajectory never carried, and
  # the only thing standing between that and a plausible wrong gradient is this
  # check: nothing about either number says which one a caller is holding.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  stand <- ladder_run(p)

  # What the run left behind is what the run introduced from and censused.
  left <- ladder_boundary_density(stand)
  patch <- ladder_as_patch(stand)
  both <- ladder_boundary_evaluations(patch, patch$ode_state, patch$ode_time)

  # Non-vacuity, and it is the whole point: if the two evaluations agreed there
  # would be no convention to get wrong and this check would pass on a model that
  # does not have the hazard.
  gap <- max(abs(both$in_uptake - both$in_field) / abs(both$in_field))
  message(sprintf("\n  the two evaluations differ by %.3e relative", gap))
  expect_gt(gap, 1e-6)

  # And the run carries the second, exactly. Tolerance is zero: this is a property
  # the implementation either has or does not.
  expect_identical(both$in_uptake, left)
  expect_false(isTRUE(all.equal(both$in_field, left, tolerance = 1e-12)))
})

test_that("the inflow condition's own parameters carry a row", {
  # The boundary node holds no ODE row, so its density adjoint has to be pushed
  # back through the condition that sets it: the adjoint at the boundary times the
  # boundary condition's own derivative. Accumulate it and never read it and the
  # accumulator is write-only -- every parameter whose only route to the census is
  # establishment then comes back short, finite and correctly signed.
  #
  # Two parameters reach the census through the inflow condition and nothing else,
  # which is what makes them the probe: the establishment constant, and the decay
  # of establishment in patch age. Neither appears in any cohort's rate equation.
  stand <- ladder_stand_introductions()
  result <- ladder_gradient_or_skip(stand)
  columns <- colnames(result$gradient)

  for (name in c("1.a_d0", "2.a_d0", "1.recruitment_decay",
                 "2.recruitment_decay")) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    peak <- max(abs(result$gradient[, name]))
    expect_gt(peak, 0)
    ladder_report_margin(paste("column", name),
                         max(abs(reference - result$gradient[, name])) / peak,
                         ladder_trajectory_agreement())
  }
})

test_that("the introduction map's whole Jacobian agrees entry by entry", {
  # Rung 5's own unit, and the one object it had no reference for. The map is the
  # pre-introduction state and the traits in, the whole widened state out; the
  # forward side seeds one tangent per input column and the reverse side is the
  # transpose the sweep runs. Both go through introduce_over, so the reference
  # traverses the forward function and not the transpose.
  #
  # Formed entirely rather than contracted, for report 08 §5.1's reason: a dot
  # product returns one number, hides an error behind a small seed component, and
  # localises to nothing when it fails. At this size the whole object fits.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29), c(0))
  stand <- ladder_run(p)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  widening <- which(diff(widths) != 0)
  times <- vapply(trajectory, function(s) s$time, numeric(1))

  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  for (k in seq_along(widening)) {
    b <- widening[[k]]
    pre <- trajectory[[b]]$state
    # The first widening introduces both species at t = 0; the second is species 1
    # alone, into a patch that already carries nodes -- which is the arrangement a
    # narrowing written as a truncation of the tail fails on.
    who <- if (k == 1L) c(1L, 2L) else 1L
    expect_equal(patch$ode_size, length(pre))

    j <- ladder_introduction_jacobian_tf24(patch, who, pre, times[[b]])
    scale <- max(abs(j$forward))
    expect_gt(scale, 0)
    gap <- max(abs(j$forward - j$reverse)) / scale
    message(sprintf("\n  widening %d: %d rows x %d columns", k, nrow(j$forward),
                    ncol(j$forward)))
    ladder_report_margin(sprintf("  introduction Jacobian, widening %d", k),
                         gap, 1e-12)
    for (i in who) patch$introduce_new_node(i)
  }
})

test_that("the state a segment resumes from is the one the introduction made", {
  # The state the segment above an introduction linearises its first step at is
  # rebuilt rather than recorded -- no record holds it -- and the sweep only
  # length-checks the rebuild. Two things must hold of it and neither was asserted:
  # the rows an introduction does not touch come through unchanged, and the rows it
  # writes are the boundary node the patch was carrying.
  #
  # Tolerance is exactly zero for both: an introduction is a push, not an
  # arithmetic operation on the rows it does not touch.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29), c(0))
  stand <- ladder_run(p)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  widening <- which(diff(widths) != 0)
  times <- vapply(trajectory, function(s) s$time, numeric(1))

  patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  for (k in seq_along(widening)) {
    b <- widening[[k]]
    pre <- trajectory[[b]]$state
    who <- if (k == 1L) c(1L, 2L) else 1L
    patch$set_ode_state(pre, times[[b]])
    invisible(patch$ode_rates)
    held <- unlist(lapply(who, function(i) patch$species[[i]]$new_node$ode_state))
    before <- vapply(seq_len(length(patch$species)),
                     function(i) patch$species[[i]]$ode_size, numeric(1))
    for (i in who) patch$introduce_new_node(i)
    post <- patch$ode_state

    stride <- length(held) / length(who)
    carried <- integer(0)
    newcomer <- integer(0)
    q <- 0L
    for (i in seq_len(length(patch$species))) {
      carried <- c(carried, q + seq_len(before[[i]]))
      q <- q + before[[i]]
      if (i %in% who) {
        newcomer <- c(newcomer, q + seq_len(stride))
        q <- q + stride
      }
    }
    carried <- c(carried, q + seq_len(length(post) - q))
    expect_identical(post[carried], pre)
    expect_identical(post[newcomer], held)
  }
})

test_that("the initial reserve channel carries a row", {
  # A germination event provisions the seedling: the initial reserve is a trait
  # times the storage capacity, so how well a seed is stocked reads the traits
  # directly. It enters nowhere else, so a dropped channel is an exact zero.
  stand <- ladder_stand_introductions()
  result <- ladder_gradient_or_skip(stand)
  at <- which(ladder_bare_traits(census_trait_names_tf24(stand)) == "a_st3")
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
  # The channel is isolated by perturbing the soil in the pre-introduction state
  # itself and introducing from each: comparing two whole runs would move the
  # answer whatever the newcomer reads, because two runs differ everywhere.
  p <- ladder_parameters(c("fast", "slow"))
  p$node_schedule_times <- list(c(0, 0.29, 0.94), c(0, 0.57))
  stand <- ladder_run(p)
  ladder_gradient_or_skip(stand)

  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  widening <- which(diff(widths) != 0)
  b <- widening[[2]]                      # a widening with a populated field
  pre <- trajectory[[b]]$state
  at <- trajectory[[b]]$time

  # The pre-introduction width is two founding nodes, one per species, which is
  # what the recorded state at this widening is a state of.
  founding_patch <- function() {
    patch <- Patch("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
    for (i in 1:2) {
      patch$set_time(0)
      patch$compute_environment()
      patch$introduce_new_node(i)
    }
    patch
  }
  newcomer_density <- function(state) {
    patch <- founding_patch()
    patch$set_ode_state(state, at)
    invisible(patch$ode_rates)
    ladder_boundary_density(patch)
  }

  # The soil block begins after every species' nodes, so its offset is read off the
  # patch the recorded state belongs to -- the founding pair, not the stand the run
  # ended at, which has three more nodes and would index past the state entirely.
  template <- founding_patch()
  expect_equal(template$ode_size, length(pre))
  moisture <- sum(vapply(seq_len(length(template$species)),
                         function(i) template$species[[i]]$ode_size,
                         numeric(1))) + 1L

  base <- newcomer_density(pre)
  # The top layer, and it has to be the top one: a plant at birth size is
  # shallow-rooted, so the deepest layer of the column does not reach it and a
  # perturbation there moves the newcomer by exactly zero. That is a property of
  # the seedling rather than of the channel, and a probe placed there would report
  # the channel missing.
  bumped <- pre
  bumped[[moisture]] <- bumped[[moisture]] * 1.05
  moved <- newcomer_density(bumped)

  gap <- max(abs(moved - base) / abs(base))
  message(sprintf("\n  a 5%% perturbation of the top layer moves the boundary node by %.3e",
                  gap))
  expect_gt(gap, 1e-9)
})

test_that("the initial condition carries no trait row on these fixtures", {
  # The last segment's state adjoint is the initial-condition term, and the sweep
  # discards the final narrowed vector. That is exact here and it is a scope
  # statement rather than a check passing: every fixture in this ladder starts
  # from bare ground, so the first segment is the empty patch and y(0) is the soil
  # alone. No trait reaches it, so there is nothing in the discarded vector.
  #
  # What the discard would lose is a run resumed from a populated exported state,
  # where founding cohorts are in y(0) and their provisioning reads a_st3. This
  # ladder has no such fixture, so the term is unrefereed rather than confirmed,
  # and saying so is the honest form: the earlier version of this test compared
  # two different runs and required the answer to move, which two different runs
  # do whatever the sweep keeps.
  stand <- ladder_stand_introductions()
  ladder_gradient_or_skip(stand)
  trajectory <- stand$store_trajectory()
  widths <- vapply(trajectory, function(s) length(s$state), numeric(1))
  # Bare ground: the first recorded width is the narrowest the run ever has, and
  # the very first step is already a widening, so nothing was integrated before the
  # first introduction. There is no cohort in y(0) for a trait to reach.
  expect_identical(widths[[1]], min(widths))
  expect_gt(widths[[2]], widths[[1]])
  message(sprintf("\n  bare ground: first width %.0f, widening at step 1",
                  widths[[1]]))
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
  columns <- ladder_bare_traits(census_trait_names_tf24(stand))
  expect_length(setdiff(declared, columns), 0L)
})

test_that("a newcomer's leaf area reaches the census through a field built without it", {
  # The fourth boundary channel, and the only one of the five whose dropped-channel
  # signature is a SHORTFALL rather than an exact zero -- so it is the one that
  # needs a reference rather than a zero test. The newcomer is placed in a field
  # built without it, so its own leaf area reaches the census through a field it
  # does not contribute to, which is a distinct path from every cohort's.
  #
  # The probe is a field-borne parameter, compared column by column against a
  # tangent of the same trajectory. On the introductions fixture every widening is
  # on the path, so this is also where the introduction boundary's own residual is
  # measured rather than inferred.
  stand <- ladder_stand_introductions()
  result <- ladder_gradient_or_skip(stand)
  columns <- colnames(result$gradient)

  # The reference has to be a derivative of the same trajectory before its
  # disagreement means anything.
  ladder_report_margin("the replay reaches this stand's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)

  # The census's own reading of the traits, so the trajectory term can be separated
  # from it. For the allometric constants the two very nearly cancel, and the
  # column's relative error is then an amplified view of the trajectory term's.
  own <- do.call(rbind, census_trait_direct_tf24(stand))
  dimnames(own) <- list(rownames(result$gradient), columns)

  for (name in c("1.a_l1", "1.a_l2", "1.k_I", "2.k_I")) {
    reference <- ladder_trajectory_tangent(
      stand, ladder_trait_direction(columns, name))$tangent
    # A shortfall, not a zero: the channel is live in both, and what is asked is
    # whether the sweep's is the same size as the reference's.
    expect_gt(max(abs(result$gradient[, name])), 0)
    r <- ladder_column_residual(result$gradient[, name], reference,
                                result$gradient[, name] - own[, name])
    borne <- ladder_bare_traits(name) %in% ladder_introduction_borne_traits()
    message(sprintf("    %-8s amplification of the trajectory term %6.2f, residual against it %.2e",
                    name, r$amplification, r$per_row_trajectory))
    ladder_report_margin(
      paste("column", name, if (borne) "(introduction-borne)" else ""),
      r$per_row,
      if (borne) ladder_introduction_residual()
      else ladder_trajectory_agreement())
  }
})

test_that("a marginal recruit sits in the stiff band, which is what makes the rest count", {
  # The fixture assertion is derived rather than chosen. The establishment
  # probability's derivative peaks at a known production, and a recruit anywhere
  # else is in a flat region where the limit below would pass on the signal being
  # small rather than on the mathematics being right.
  #
  # At the shipped establishment constant this fixture's recruits sit at 36 and 21
  # times that production, with an establishment probability of 0.998 -- a
  # comfortable seedling, and the wrong place to ask the question.
  stand <- ladder_stand_marginal_recruit()
  ladder_gradient_or_skip(stand)
  stiffness <- ladder_recruit_stiffness(stand)
  establishment <- ladder_recruit_establishment(stand)
  message(sprintf("\n  recruit at %.3f of the peak production, establishment %.5f",
                  stiffness, establishment))
  expect_true(all(stiffness > 1 / 3 & stiffness < 3))
})

test_that("a marginal recruit's gradient exists, is finite, and tends to zero", {
  # Both the cumulative mortality and the log density diverge logarithmically as
  # establishment goes to zero, and the sensitivity of the log density diverges
  # with them. But the log density reaches everything downstream only through the
  # stem number, whose own sensitivity tends to zero, and the census carries stem
  # number linearly. So the finite answer is reached by seeding the quantity that
  # has one, and a recruit that cannot pay for itself contributes at second order
  # -- exactly zero is the correct first-order answer.
  #
  # The probe is the establishment decay, which reaches the census through the
  # boundary condition and nothing else, so its column is the recruit's own
  # contribution rather than a stand-wide one.
  scales <- c(30, 60, 120, 400, 2000)
  establishment <- numeric(length(scales))
  decay <- numeric(length(scales))
  for (k in seq_along(scales)) {
    stand <- ladder_stand_marginal_recruit(scales[[k]])
    result <- ladder_gradient_or_skip(stand)
    expect_true(all(is.finite(result$gradient)))
    establishment[[k]] <- ladder_recruit_establishment(stand)[[1]]
    columns <- colnames(result$gradient)
    decay[[k]] <- max(abs(result$gradient[
      , ladder_bare_traits(columns) == "recruitment_decay"]))
    message(sprintf("  establishment %.6f   |recruitment_decay| %.4e",
                    establishment[[k]], decay[[k]]))
  }

  # Non-vacuity: the sweep has to actually approach the limit, or "tends to zero"
  # is a statement about one point.
  expect_lt(min(establishment), 1e-3)
  expect_gt(max(establishment) / min(establishment), 1e3)

  # It tends to zero, and it does so from above at every step rather than only
  # between the endpoints.
  expect_true(all(diff(decay) < 0))
  expect_lt(decay[[length(decay)]] / decay[[1]], 1e-2)
})
