# Dead code found while tracing the gradient

Noted rather than removed. Each entry says what makes it dead and what would have
to be true for it to be live again, so the removal is a decision rather than a
rediscovery.

## Removed

- **`inst/include/plant/collar_census.h`** — a census of collar-solve termination
  classes, gated behind `PLANT_COLLAR_CENSUS`. No call sites anywhere in the
  workspace: the leaf moved to phylloptim and the census stayed behind. It is
  superseded by `phylloptim::Leaf::operating_point_kind()`, which classifies by
  the branch taken, per solve, reset at the top — the same thing this header's
  own comment describes, done properly. Its comment is also the origin of the
  "a loose `GSS_tol_abs` reports nearly every solve as pinned" warning that the
  `GSS_tol_abs` known issue cites.

## Found, not removed

- **`inst/include/plant/optimize.h`** — Brent's method and a golden-section
  maximiser. Dead on three counts: nothing includes it; `brent_fmin` and
  `golden_section_max` are referenced nowhere in the tree; and its own comment
  says it exists for "the leaf hydraulic solvers", which both moved to phylloptim
  and retired golden-section in favour of solving the first-order condition.
  Live again only if something in plant needs a 1-D minimiser of its own.

- **The `!birth_date` branches inside the reverse pass.** The gradient is scoped
  to the birth-date coordinate and refuses the height one at entry, so the
  weight-derivative terms in `Species::compute_competition_and_slope_adjoint` —
  the `out[upper].height += edge; out[k].height -= edge;` pair and the guard
  above them — cannot be reached on any differentiated path. They are the
  transpose of a coordinate the sweep will not run. The forward model still
  supports both coordinates, so only the *adjoint* copies are dead. Live again
  if the gradient is ever extended to the height coordinate, which report 05
  §6.1 and report 06 §11 both argue against.

  *Not verified:* I did not confirm where the height-coordinate refusal is
  enforced. Confirm that before removing anything here — if the refusal is at
  one entry point rather than all of them, these branches are reachable.

- **`boundary_node_adjoints` has three write-only members.** `.area_leaf`,
  `.height` and `.extinction` are written from four sites in `patch.h` and read
  from none; only `.density_in_field` and `.density_in_uptake` are consumed, by
  `boundary_condition_adjoint`. Each is accounted for — `.height` is the imposed
  `dh0/dphi = 0`, and `.area_leaf`/`.extinction` are already consumed by
  `light_reduction_trait_adjoint`, which runs for the boundary node before the
  branch that writes them — so these are redundant writes rather than a dropped
  channel. Worth deleting because a write-only accumulator is indistinguishable
  by inspection from a channel that was forgotten.

- **`Node::set_log_density_rate`** (`node.h:121`) — declared and called from
  nowhere in the tree. The density rate is assigned by another route.

- **`Species::growth_rate_gradient`** — the compression term, `-dg/dh` differenced
  across neighbouring nodes. Its only caller is inside
  `if (internals::transport_census_active())`, a diagnostic behind an environment
  variable, so it is off the rate path entirely. That is correct on the
  birth-date coordinate, where the compression term does not exist; it is
  reachable again only through that diagnostic or through the height coordinate.
  Checked because a stencil differencing an active quantity across neighbours
  would have matched the a_l1/a_l2 defect's signature exactly — per species,
  switching on at the second cohort. It does not, because it never runs.

## Live, despite looking otherwise

- **`inst/include/plant/leaf_model.h`** — a 63-line compatibility shim aliasing
  `plant::Leaf = phylloptim::Leaf` so the generated RcppR6 and RcppExports glue,
  which names `plant::Leaf` throughout, keeps compiling. Deliberate and
  documented. It does pull in constants that nothing reads, on the stated grounds
  that they stay reachable as `plant::<name>` "for anything that reaches for them
  later"; only `kg_per_mol_h2o` is used today.
