# Run the checks on the reverse-mode gradient, and print where they stand.
#
#   Rscript scripts/run-gradient-ladder.R                 every file
#   Rscript scripts/run-gradient-ladder.R one-cohort       one file
#
# Run from the package root, and build optimised first: an unoptimised build
# makes the trajectory fixtures several times slower.
#
#   make
#
# ---------------------------------------------------------------------------
#
# WHAT THESE CHECKS ARE, AND WHY THERE ARE SO MANY OF THEM
#
# The claim under test is that the reverse sweep is the exact transpose of the
# forward run. That claim cannot be checked by a finite difference, because
# differencing is the thing the sweep replaces -- and near a coincidence a
# difference does not merely lose accuracy, it fails to converge. So the claim is
# checked against references that share no code with the sweep, and the files
# below are grouped by which reference each one uses.
#
# There are four references, and none of them alone is enough:
#
#   1. A FORWARD-MODE TANGENT of the same forward source. Exact -- no step size,
#      no truncation -- and it traverses the forward reductions while touching
#      none of the transposes under test. One seed gives one exact column.
#      ⚠️ It cannot see a wrong SUPPLIED row: the tangent runs the same supplied
#      numbers through the same attachment, so a bad row makes both routes wrong
#      identically. That is what reference 4 is for.
#
#   2. THE WHOLE JACOBIAN, formed entry by entry, at one cohort. A contraction
#      against a seed returns one number, hides an error behind a small seed
#      component, and localises to nothing when it fails. Affordable only at the
#      smallest fixture, which is why that fixture exists.
#
#   3. A DIFFERENCE OF WHOLE RUNS in plain double. The only reference that shares
#      no arithmetic at all with the sweep, and therefore the only one that can
#      catch a common-mode error. Coarse, so it bounds rather than pins.
#
#   4. THE FORWARD MODEL REBUILT FROM ITS PARAMETERS, which is the one place
#      references 1 and 2 are blind together. A rebuild runs preparation, so the
#      leaf is CONSTRUCTED with the moved trait rather than handed it afterwards,
#      and differencing the rates then traverses the leaf's forward solve while
#      touching none of the derivative code the supplied rows come from.
#      ⚠️ It is needed because differencing the recorded step cannot referee a
#      supplied row at all: the recorded expression is the value plus a sum of
#      partials times brackets that are each exactly zero, so the difference comes
#      back EXACTLY zero on precisely the columns a supplied row occupies --
#      whether the row is right, wrong, or absent.
#
# The fixtures form a progression -- the "ladder" the filenames name -- from the
# smallest stand whose Jacobian can be formed whole, through accumulation across
# cohorts and species, to a stand whose state vector changes width. Each step
# adds one mechanism and uses the strongest reference its size still permits.
#
# ⚠️ A MARGIN IS NOT A PASS. Every check reports how much of its tolerance budget
# it used, because a check passing at three-quarters of budget is a check about to
# stop working. And a suite of margins says nothing about whether the checks would
# notice a defect at all, which is what the fault injections establish: two of the
# first three injections tried here failed to fail.
#
# ⚠️ A SKIP IS NOT A PASS EITHER. A skip is these checks naming something they
# cannot yet ask. Read the skip count beside the failures every time.

suppressMessages({
  library(odelia)
  pkgload::load_all(".", quiet = TRUE)
  library(testthat)
})

# What each file claims, in the order a reader should meet them. The key is the
# filename stem after `test-gradient-ladder-`.
claims <- c(
  # No reference needed: these hold or the model is wrong.
  "model-invariants"    = "the model's own invariants, before any reference exists",
  "identity"            = "the sweep gives the same answer however it is decomposed",

  # Reference 2, the whole Jacobian, at the only size that permits it.
  "one-cohort"          = "one cohort, and the Jacobian formed entirely",
  "switches"            = "one switch per route a parameter takes to a census",

  # Reference 1, the forward tangent, over trajectories.
  "two-species"         = "two species, two cohorts each, over a trajectory",
  "columns"             = "each shortlisted trait's own column, against a tangent",
  "introductions"       = "introductions, where the state vector changes width",
  "first-range"         = "the range of the recording below the first introduction",
  "recruit"             = "the inflow boundary, and how often its own row enters",

  # Reference 4, a rebuilt forward model, on the one seam the others cannot see.
  "factorisation"       = "the one seam: the rows the reference and the sweep share",

  # Reference 3, a difference of whole runs.
  "whole-run-difference" = "the sweep against a difference of whole runs, five regimes",
  "declared-zero"       = "pricing a zero the model declares rather than computes",

  # What the checks do when the model declines to answer, and whether they bite.
  "sweep"               = "whether the sweep runs at all, and on which channels",
  "injection"           = "whether these checks bite, established by breaking what they watch",
  "production-scale"    = "what a production-length stand violates, recorded not enforced")

# ⚠️ Four more gradient files sit outside this progression and outside this
# script, because they check the surface rather than the sweep:
# `test-finite-difference.R` (the helper this file's name used to collide with),
# `test-gradient-demo.R`
# (the staging study's own helpers), `test-gradient-incidence.R` (how often each
# regime is met) and `test-gradient-parity.R` (that every state the forward model
# answers for, the reverse either answers or names). Run those with the ordinary
# test runner.

dir <- file.path("tests", "testthat")

# ⚠️ MAINTAINED, AND CHECKED, because a list read off a directory goes silently
# incomplete when the directory grows. A file with no entry here would otherwise
# never run under this script and would look like a suite that passes.
on_disk <- sub("^test-gradient-ladder-", "",
               sub("\\.R$", "",
                   basename(Sys.glob(file.path(dir, "test-gradient-ladder-*.R")))))
undescribed <- setdiff(on_disk, names(claims))
if (length(undescribed) > 0L) {
  stop("These files have no entry in `claims`, so this script would not run ",
       "them: ", paste(undescribed, collapse = ", "),
       ". Add one saying what the file claims.")
}
missing <- setdiff(names(claims), on_disk)
if (length(missing) > 0L) {
  stop("These entries in `claims` name no file: ",
       paste(missing, collapse = ", "))
}

selected <- commandArgs(trailingOnly = TRUE)
if (length(selected) == 0L) selected <- names(claims)
unknown <- setdiff(selected, names(claims))
if (length(unknown) > 0L) {
  stop("No such file: ", paste(unknown, collapse = ", "),
       ". Available: ", paste(names(claims), collapse = ", "))
}

results <- list()
for (name in selected) {
  path <- file.path(dir, paste0("test-gradient-ladder-", name, ".R"))
  cat(sprintf("\n===== %s -- %s\n", name, claims[[name]]))
  out <- testthat::test_file(path, package = "plant", reporter = "summary")
  frame <- as.data.frame(out)
  results[[name]] <- c(
    pass = sum(frame$passed),
    fail = sum(frame$failed),
    skip = sum(frame$skipped),
    error = sum(frame$error))
}

cat("\n\n===== where the checks stand =====\n")
table <- do.call(rbind, results)
print(table)

blocked <- rownames(table)[table[, "skip"] > 0]
if (length(blocked) > 0) {
  cat("\nSkips name something these checks cannot yet ask, and they are not\n",
      "passes. Skipped: ", paste(blocked, collapse = ", "), "\n", sep = "")
}
