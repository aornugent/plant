# AD emergent-gradient roadmap (#472 scope B / #537)

Status as of the `spike-ff16-scm-emergent` work: the live-SCM two-pass emergent
trait-gradient framework is built and validated end-to-end for FF16 — both shading
models (crown-top, deep-crown), the real `offspring_production`, the full 28-trait
sweep in one reverse pass, the differentiated establishment filter, and both the
frozen-resident (invasion) and active-knot (resident-reshaping) light treatments.
Every piece is checked against finite differences and faithful to the live SCM
(cohort heights to ~1e-14, `offspring_production` to ~1e-13). See `scripts/ad_*.R`.

What it is NOT yet: landed on `develop`, CI-tested end-to-end, or exposed as a
callable calibration API. This roadmap is the plan from here.

## Phase A — Land the foundation (gate; not coding)

- Get **PR #541** (`spike-ff16-hierarchy`, the `<T,E,S>` hierarchy + Node exact
  gradient) reviewed and merged to `develop`.
- Rebase: `git rebase --onto develop spike-ff16-hierarchy spike-ff16-scm-emergent`,
  then reinstall against the clean develop baseline.
- **Plan:** keep building C + D on the spike stack while A is in review; the reviewer
  looks at A (#541) and the D refinement PR together.

## Phase C — Production calibration entry point (in progress)

Turn the validated machinery into something a calibration loop calls with no
on-the-fly compilation.

- Move the two-pass reverse-mode replay from sourceCpp into a compiled `.cpp` in
  `src/` and expose it via RcppR6 (the way `Individual$growth_rate_gradient_exact`
  / `growth_rate_gradient_height_ad` are compiled into `plant.so`).
- Entry point: a gradient of an emergent output (`offspring_production`) w.r.t. a set
  of traits, returned as a named vector — usable directly as a calibration objective
  gradient (many traits → one objective → one reverse sweep).
- This makes the headline result CI-testable in plain R (no odelia-link / BH-less
  skip), closing the main pre-merge gap.

## Phase D — Fidelity refinements (in progress; PR reviewed with A)

- Per-**RK-stage** stand harvest → bit-exact active resident light (currently
  per-step ⇒ ~1e-4 reconstruction).
- Multi-species stand harvest (currently species-0 only).
- `d(height_0)/d(trait)` → completes establishment differentiation for seedling-size
  traits (`omega`, `lma`; `height_0` is derived in `prepare_strategy`).

## Phase E — Scientific payoff (the point; eco-evo priority)

- **Selection gradients for adaptive dynamics:** the frozen-resident invasion-fitness
  gradient IS the selection gradient — gradient ascent in trait space / locating
  evolutionary singular strategies. Needs no new machinery.
- **Gradient-based calibration:** fit FF16 traits to data through the emergent
  gradients (the 28-trait sweep against a likelihood).
- **Acclimation** (#406/#537/#527): the unified whole-plant growth objective gradient.

## Phase F — Full AD for TF24 (new)

TF24 (hydraulics + flexible allometry + NSC storage) is the harder strategy: its RHS
carries a **root-find** (leaf `psi_stem → ci`) and **splines** (vulnerability curve,
environment), where FF16 has neither. The groundwork exists:
- A2 / **#539 (MERGED)**: the leaf-level gradient — d(profit)/d(root-collar psi) via
  forward-mode AD + the implicit function theorem at the `psi_stem→ci` root-find.
- **odelia #32**: scalar-templated differentiable spline.
- The live-SCM two-pass replay machinery (harvest bindings, generic Cash-Karp stepper)
  is strategy-agnostic — what is TF24-specific is the rate kernels.

Mirror the FF16 path, simplest-config-first (well-watered before drought):
- **F1 — leaf/net-production kernel:** scalar-template the TF24 assimilation + net
  production. Carry the gradient through the `psi_stem→ci` root-find via IFT (the A2
  precedent: seed the derivative at the converged root, do not differentiate the
  solver iterations). Delegate the double methods to it (faithfulness vs the TF24
  reference tests).
- **F2 — demographic rate fill + NSC:** growth / fecundity / mortality, plus the NSC
  storage dynamics and the (drought) soil-water coupling. Larger/different state
  vector than FF16's 5.
- **F3 — emergent gradient:** reuse the generic two-pass replay (the harvest bindings
  already instantiate for TF24/TF24_Env). Validate replay vs live TF24 SCM heights and
  d(emergent)/d(trait) vs FD, well-watered first, then with drought.

## Sequencing

A is the gate but does not block C/D (built additively on the spike stack; A + the D
PR reviewed together). **C before D** — get a callable API in front of the validated
machinery before chasing fidelity; current fidelity already suffices for the Phase-E
selection gradient. **E is the goal** — let a concrete selection-gradient or
calibration use case drive exactly which D refinements are worth doing. **F (TF24)**
is the second big strategy; it reuses the whole framework and is gated only by its own
kernels, so it can proceed in parallel once C gives the API pattern to follow.
