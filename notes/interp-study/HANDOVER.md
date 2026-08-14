# Picking this up: the light field's knot placement

What is landed, what is open, and what has to be rebuilt because it lived in a
scratch directory that does not survive the session.

---

## 1. Read this first, in this order

1. [`../light-field-knot-placement.md`](../light-field-knot-placement.md) — the
   design and every measurement behind it. §5 is what was built, §6 is what was
   refuted and why, **§8.5 is the `develop` comparison and the frontier
   argument**, §9 is what would falsify it.
2. `NEWS.md`, top two entries — the R-facing change and the migration.
3. `git log --oneline ec8e59f1..ad/light-field-fixed-grid` — six commits, each
   message states what it changed and why.

The corpus reports this rests on: **03 §3.1** (why the interpolant is Hermite),
**03 §3.3** (the passive-position treatment this is about), **01 §7 item 8** (why
a state-dependent knot count is ruled out), **05 §5** (the recorded block's
inputs), **08 §4.8** (two paths agreeing because they share a multiplier — the
mechanism that hid the defect).

---

## 2. Where the work is

| what | where | state |
|---|---|---|
| the change | `plant` branch `ad/light-field-fixed-grid` | **2 commits unpushed** at time of writing |
| its worktree | `plant/.claude/worktrees/fixed-grid` | built at `-O2` |
| odelia's test | `odelia` branch `ad/hermite-multiplied-grid` | pushed, 1 commit |
| the measurement scaffolding | `plant` branch `ad/interp-design` | **unpushed**, worktree `plant/.claude/worktrees/interp-design` |
| `develop`, for comparison | worktree `plant/.claude/worktrees/develop` | built, needs odelia v0.2.1 (§4) |

`ad/interp-design` is worth keeping: it carries `ResourceGridPolicy`, which
selects the knot placement **at run time**, so one binary can time and referee
every candidate. Every placement comparison in the design note was made with it.
It is scaffolding, not a shipped control, and must not be merged.

**Check before anything else:**

```sh
cd plant && git log --oneline -6 ad/light-field-fixed-grid
git status --short                      # expect clean
git log --oneline origin/ad/light-field-fixed-grid..ad/light-field-fixed-grid
```

---

## 3. The base has moved — rebase, and regenerate rather than merge

The branch was cut from `ec8e59f1`. `ad/v3-forward` has since advanced to
`7015a1c9` ("Re-derive the seed's dependent aux where its height is written"),
which was the staged WIP when this work started.

```sh
cd plant/.claude/worktrees/fixed-grid
git rebase ad/v3-forward
```

**`R/RcppExports.R` and `src/RcppExports.cpp` will conflict, and they are
generated.** Both sides edit them: this branch because `ResourceSpline`'s
constructor changed, `7015a1c9` because it adds exports. Do not hand-merge —
take either side and regenerate:

```sh
Rscript -e 'library(methods); RcppR6::RcppR6()'   # also runs compileAttributes
```

`inst/RcppR6_classes.yml` is the real source for the `ResourceSpline`
constructor and is hand-edited; that one **is** a genuine conflict if both sides
touched it (they should not have).

---

## 4. What must be rebuilt (it lived in scratch)

**The odelia v0.2.1 private library, needed to build `develop`.** `develop` pins
`odelia@v0.2.1`, and the installed odelia has since removed
`odelia::ode::iterator`, so `develop` will not compile against it. Rebuild into a
private library so the installed one is left alone:

```sh
cd odelia && git worktree add .claude/worktrees/v021 v0.2.1     # if gone
mkdir -p /tmp/rlib-dev
Rscript -e 'install.packages("odelia/.claude/worktrees/v021", repos=NULL,
                             type="source", lib="/tmp/rlib-dev")'

cd plant/.claude/worktrees/develop && rm -f src/*.o src/*.so
R_LIBS_USER=/tmp/rlib-dev:$HOME/R/x86_64-pc-linux-gnu-library/4.6 \
R_MAKEVARS_USER=$PWD/../fixed-grid/notes/interp-study/Makevars-O2 \
Rscript -e 'pkgbuild::compile_dll(".", compile_attributes = FALSE, debug = FALSE)'
```

Everything else that mattered is now in this directory: `Makevars-O2`, the
harness, and the test runners. The scripts resolve their siblings relative to
themselves, so they run from the repo root.

**Build at `-O2`.** `pkgbuild::compile_dll` defaults to `-O0` and appends it
after any user flags, so a timing taken without this measures the debug build.

```sh
R_MAKEVARS_USER=$PWD/notes/interp-study/Makevars-O2 \
Rscript -e 'pkgbuild::compile_dll(".", compile_attributes = FALSE, debug = FALSE)'
```

---

## 5. Re-running the measurements

`lib-field.R` and `d1-lib.R` are the two harnesses. Everything referees against
the model's own reduction (`patch$compute_competition_and_slope`), so no state
crosses between builds.

| question | script | notes |
|---|---|---|
| does the harness reproduce the shipped interpolant? | `h3-validate.R` | **run this first**; expect 1.1e-16 |
| accuracy and dropped channel, every placement | `h4-placements.R` | |
| dropped channel per cohort; append inertness; canopy sweep | `h5-sweep.R` | |
| end to end: tangent vs a difference of the rates | `h7-endtoend.R` | the 2127x -> 1.2x result |
| is it the count or the break crossing? | `h9-crossing.R` | |
| what a placement costs | `h6-cost.R`, `d4-cost.R`, `d5-where.R` | |
| convergence, so a re-bless can be justified | `h14-converge.R`, `h21` (in log) | |
| **through stand development, per build** | `d2-early.R` | `PLANT_BUILD=<worktree> LABEL=... ` |
| **thinning, all placements in one binary** | `d3-thinning.R` | the table that changed the spacing |
| the dyadic alternative | `d6-dyadic.R`, `d7-level.R` | §7 below |

Test runners: `run-all.R` (whole suite, ~40 min), `run-subset.R` (`FILES=a,b`),
`show-fail.R` (prints failure detail). **Set `NOT_CRAN=true`** or the
`model-version` drift guard silently skips — it is `skip_on_cran()`, and it is
the one test that catches an unbumped `scientific_version`.

---

## 6. State of play

**Landed and verified.** Knots at `k * spacing`, spacing 0.05 m for FF16/TF24 and
0.025 m for K93 (≈ the canopy a run reaches / 350). The dropped position channel
is exactly zero; the tallest cohort's height column went from 2127x the other
height columns to 1.2x. Full suite **3625 pass / 5 fail** at the previous spacing,
with all 5 pre-existing (`test-mutant.R` 2, `test-stochastic-patch.R` 3 — they
fail identically on the shipped placement, verified).

**Unfinished at session end:** the confirmatory full suite at the corrected
spacing (0.05) had not completed. The affected files were re-run clean
individually — canopy-methods 282, ff16 55, k93 21, environment 33, scm 140,
patch 167, ladder rung3 41 — but **the whole-suite number at 0.05 is not in
hand.** Run `run-all.R` first thing.

**Re-blessed:** `deep-crown reproduces the baseline SCM result` 16.8846 ->
16.8954 (limit 16.8961), `offspring arrival` in `test-strategy-ff16.R`.
`scientific_version` bumped FF16 2->3, K93 2->3, TF24 9->10, snapshot accepted.

---

## 7. Open, in the order they matter

1. **The frontier question is not closed.** A lattice does not dominate
   canopy-tied knots at matched cost — it wins late and loses early. What it buys
   is that the adjoint differentiates the function the forward model computes.
   §8.5 of the design note has the numbers and the argument.

2. **The dyadic lattice is the measured route to recovering the cost**, and it is
   not built. Every `2^m`-th point of one fixed lattice, level chosen from the
   canopy: position channel still exactly zero, count bounded at 182, beats
   canopy-tied knots at every canopy above 1.5 m, ~1.03x run cost against the
   0.05 lattice's 1.55x. **Its cost is four field discontinuities of 1.5e-05 to
   3.8e-04 of range per run**, where a uniform lattice has only a kink. Removing
   that jump — a blend across the level change — is the work. `d6-dyadic.R` and
   `d7-level.R` are the measurements to build against.

3. **The shaded-stand non-finite tangent bounds what can be refereed.** Clean at
   a minimum read light of 0.203, NaN at 0.143 and below; 38 of 73 rate rows
   wholly NaN, the same count at three shade levels, so it is a switch. Placement
   independent — present on both. It means the design is **unrefereed end to end
   in deep shade**, which is where the field matters most.

4. **Report 01 §3's "peak is flat in run length" is now false** via the canopy:
   one cohort's recording is linear in the knot count at 480 slots per knot. Flat
   in *cohort* count still holds exactly (identical at 1, 2, 4, 8). That is a
   corpus claim and amending it is the corpus owner's call.

5. **Two placement-independent column disagreements** sit beside the ladder: a
   `storage` column at 2.9e-01 and two `log_density` columns at ~2.5e-03,
   identical under every placement. Most likely report 08 §4.7's trap two — the
   difference's step is absolute for a small state — not a new defect, but
   unconfirmed.

---

## 8. Things that will bite

- **`pkill -f run-all.R` kills the run you just started in the same command.**
  Use `nohup ... & disown` and kill by PID.
- **The `model-version` drift guard is `skip_on_cran`.** Without `NOT_CRAN=true`
  a `scientific_version` that should have been bumped passes silently.
- **`patch$environment` returns a copy**, so mutating it in R does not reach the
  patch. A test that sets a fixed environment through it is testing a temporary.
- **`max_environment_height()` is not the canopy top** any more — the field
  reaches one knot past the tallest cohort. It is a bound on the field.
- **The interp-design build crashes with `malloc(): invalid size`** in fixed mode
  under the gradient ladder: it predates the `knot_count()` fix. Use the
  fixed-grid worktree for anything touching the recorded block.
