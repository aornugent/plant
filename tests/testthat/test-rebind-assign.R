# A patch reaches another scalar two ways: rebind_from builds one, assign_from
# writes into one that exists. Every comment in both packages says they leave the
# patch holding the same thing, and until this file nothing checked it -- so when
# they drifted, what said so was the gradient ladder, twenty minutes away, in
# columns that had no obvious connection to a copy.
#
# The two do not fail the same way. A rebind returns a fresh object, so a member
# it does not write is default-constructed. An assignment writes into one that
# exists, so a member it does not write keeps what it had. Every member is
# therefore a place they can differ, and the check has to be over the whole
# object rather than over the ones someone remembered.

test_that("a rebound patch and an assigned one hold the same thing", {
  p <- ladder_parameters(c("fast", "slow"))

  # Somewhere with cohorts and a built field, not bare ground: an empty patch
  # agrees trivially, because the members that drift are the ones a run fills.
  scm <- SCM("TF24", "TF24_Env")(p, Environment("TF24"), ladder_control())
  scm$run()
  patch <- scm$patch

  got <- ladder_rebind_matches_assign_tf24(patch)

  # Widths first. A width that disagrees is the failure that reaches furthest --
  # n_cohort_reads is 2 * knot_count + soil depths, so it sets the width of every
  # recorded block, and a difference here misshapes the whole reverse pass rather
  # than moving a number in it.
  expect_identical(got$ode_size[[1]], got$ode_size[[2]])
  expect_identical(got$aux_size[[1]], got$aux_size[[2]])
  expect_identical(got$n_cohort_reads[[1]], got$n_cohort_reads[[2]])
  expect_identical(got$n_parameters[[1]], got$n_parameters[[2]])

  # Then values, and exactly: both sides carry the same doubles through the same
  # conversion, so anything above zero is a member one path writes and the other
  # does not.
  expect_identical(got$state_gap, 0)
  expect_identical(got$parameter_gap, 0)
})
