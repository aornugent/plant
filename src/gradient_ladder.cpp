#include <limits>
#include <plant.h>

// The ladder's references for the census gradient, and the instrumentation
// that reports how a run was cut and what regimes it met.

// The references the gradient is checked against, and the two objects they are
// checked at. The two objects referee different claims, and reading them as one
// overstates both.
//
// The main reference is a forward-mode tangent of the same forward source. It is
// exact -- no step size, no truncation -- and it traverses the forward
// reductions while none of the transposes under test are on its path. One seed
// gives one exact Jacobian column, which is what makes forming a whole Jacobian
// affordable at one cohort.
//
// One right-hand-side evaluation relates the whole ODE state to its own rates
// through both reductions. Taking it forward and backward checks the sweep, and
// there the reference and the object under test share no code.
//
// The recorded cohort block relates one individual's states, the field's knot
// values and slopes, the soil potentials and the traits to its rates, its density
// rate and its per-layer draws. Nothing in production records it; what it
// referees is the recording. Taking it forward and backward checks that both modes
// traverse the same recording, and it cannot check the number a recorded row
// carries: both modes read that number from the same place, so they agree
// whether it is right or wrong. The plain-double difference below is the only
// instrument here that can see a wrong supplied row, because it runs the leaf's
// own solve instead of reading the row.

namespace {

using strategy_type = plant::TF24_Strategy<double>;
using environment_type = plant::TF24_Environment<double>;
using patch_type = plant::Patch<strategy_type, environment_type>;
using patch_handle = plant::RcppR6::RcppR6<patch_type>;

// Both scalars come from the one place each is named -- patch.h for the tangent,
// odelia for the adjoint -- rather than being restated here, and the two things
// done to a tangent come from the same place.
using plant::tangent;
using plant::seed_direction;
using plant::derivative_along;
using tangent_strategy = plant::at_scalar<strategy_type, tangent>;
using tangent_environment = plant::at_scalar<environment_type, tangent>;
using tangent_individual = plant::Individual<tangent_strategy,
                                             tangent_environment>;

using adjoint = odelia::ode::active_scalar<double>;
using adjoint_patch = plant::Patch<plant::at_scalar<strategy_type, adjoint>,
                                   plant::at_scalar<environment_type, adjoint>>;

// One right-hand-side transpose, taken by the call the sweep takes: odelia rebinds
// the patch to the adjoint scalar, records derivs() and sweeps the recording per
// seed. Returns the recording's size. The tape is the caller's, as it is the
// solver's in a sweep, so a caller repeating the call does not pay for it twice.
size_t rhs_adjoint(patch_type& patch,
                   odelia::ode::adjoint_tape<double>& tape,
                   const std::vector<double>& lambda_dydt,
                   std::vector<double>& lambda_state,
                   std::vector<double>& trait_adjoint) {
  const size_t n = patch.ode_size();
  plant::util::check_length(lambda_dydt.size(), n);
  std::vector<double> state(n);
  patch.ode_state(state.begin());

  const odelia::ode::adjoint_rows seeds =
    odelia::ode::adjoint_rows::one_row(lambda_dydt);
  odelia::ode::adjoint_rows swept;
  odelia::ode::adjoint_rows rows(1, patch.trait_adjoint_size());
  odelia::ode::active_system<std::decay_t<decltype(patch)>> active{patch, tape};
  const size_t recording = odelia::ode::rates_adjoint(
    active, state, patch.time(), seeds, swept, rows);
  lambda_state.assign(swept[0].begin(), swept[0].end());
  trait_adjoint.assign(rows[0].begin(), rows[0].end());
  return recording;
}

// Which species and which of its nodes a flat node index names, in the order
// the block loop and the state vector both visit them.
struct node_address { size_t species; size_t node; };

node_address locate(const patch_type& patch, size_t flat) {
  size_t seen = 0;
  for (size_t i = 0; i < patch.size(); ++i) {
    const size_t n = patch.at_species(i).size();
    if (flat < seen + n) {
      return node_address{i, flat - seen};
    }
    seen += n;
  }
  plant::util::stop("node index past the end of the patch");
  return node_address{0, 0};
}

// The block as the sweep records it, standing ready to be evaluated at a seeded
// input vector. The strategy and environment are rebound once and copied per
// evaluation, because a rebind runs the strategy's preparation and preparing at
// an active scalar is refused.
struct tangent_block {
  tangent_strategy strategy_template;
  tangent_environment environment_template;
  std::vector<double> in;
  size_t n_out;

  tangent_block(const patch_type& patch, const node_address& at)
    : strategy_template(
        patch.at_species(at.species).strategy_ptr()->rebind_from<tangent>()),
      environment_template(patch.r_environment().rebind_from<tangent>()) {
    // The knot positions come across because a rebind carries the field's data
    // and not its grid.
    environment_template.light_availability.spline.set_nodes(
      patch.r_environment().light_availability.spline.knots());

    typename tangent_strategy::ptr strategy =
      std::make_shared<tangent_strategy>(strategy_template);
    tangent_individual probe(strategy);
    in.resize(probe.block_input_size(environment_template));
    n_out = probe.block_output_size(environment_template);
    patch.at_species(at.species).node_at(at.node).individual
      .block_inputs(in.begin(), patch.r_environment());
  }

  // One evaluation with column `seed` carrying a unit tangent. Returns the
  // outputs' tangents, which is one exact Jacobian column.
  std::vector<double> column(int seed) {
    typename tangent_strategy::ptr strategy =
      std::make_shared<tangent_strategy>(strategy_template);
    tangent_individual individual(strategy);
    tangent_environment environment = environment_template;

    std::vector<tangent> x(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
      x[i] = in[i];
    }
    if (seed >= 0) {
      seed_direction(x[static_cast<size_t>(seed)], 1.0);
    }

    individual.set_block_inputs(x.begin(), environment);
    individual.compute_rates(environment);
    std::vector<tangent> y(n_out);
    individual.block_outputs(y.begin(), environment);

    std::vector<double> out(n_out);
    for (size_t j = 0; j < n_out; ++j) {
      out[j] = derivative_along(y[j]);
    }
    return out;
  }

  std::vector<double> value() {
    typename tangent_strategy::ptr strategy =
      std::make_shared<tangent_strategy>(strategy_template);
    tangent_individual individual(strategy);
    tangent_environment environment = environment_template;
    std::vector<tangent> x(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
      x[i] = in[i];
    }
    individual.set_block_inputs(x.begin(), environment);
    individual.compute_rates(environment);
    std::vector<tangent> y(n_out);
    individual.block_outputs(y.begin(), environment);
    std::vector<double> out(n_out);
    for (size_t j = 0; j < n_out; ++j) {
      out[j] = odelia::util::to_passive(y[j]);
    }
    return out;
  }
};

Rcpp::NumericMatrix to_matrix(const std::vector<std::vector<double>>& rows,
                              size_t n_col) {
  Rcpp::NumericMatrix out(static_cast<int>(rows.size()),
                          static_cast<int>(n_col));
  for (size_t i = 0; i < rows.size(); ++i) {
    plant::util::check_length(rows[i].size(), n_col);
    for (size_t j = 0; j < n_col; ++j) {
      out(static_cast<int>(i), static_cast<int>(j)) = rows[i][j];
    }
  }
  return out;
}

} // namespace

// The block's inputs, in the order block_inputs writes them: the individual's
// own states, the field's knot values then its knot slopes, the soil potentials,
// and the strategy's differentiable parameters.
// [[Rcpp::export]]
std::vector<std::string> ladder_block_input_names_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                       int node) {
  const patch_type& patch = *obj_;
  const node_address at = locate(patch, static_cast<size_t>(node - 1));
  std::vector<std::string> out = strategy_type::state_names();
  const size_t n_knot =
    patch.r_environment().light_availability.spline.knots().size();
  for (size_t i = 0; i < n_knot; ++i) {
    out.push_back("light_value_" + plant::util::to_string(static_cast<int>(i + 1)));
  }
  for (size_t i = 0; i < n_knot; ++i) {
    out.push_back("light_slope_" + plant::util::to_string(static_cast<int>(i + 1)));
  }
  const int n_layer = patch.r_environment().get_soil_number_of_depths();
  for (int i = 0; i < n_layer; ++i) {
    out.push_back("psi_soil_" + plant::util::to_string(i + 1));
  }
  for (const std::string& n :
       patch.at_species(at.species).strategy_ptr()->ad_parameter_names()) {
    out.push_back(n);
  }
  return out;
}

// The block's outputs: the strategy's rates, the density rate, then one
// consumption rate per resource.
// [[Rcpp::export]]
std::vector<std::string> ladder_block_output_names_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  std::vector<std::string> out;
  for (const std::string& n : strategy_type::state_names()) {
    out.push_back(n + "_dt");
  }
  out.push_back("log_density_dt");
  const int n_layer = patch.r_environment().get_soil_number_of_depths();
  for (int i = 0; i < n_layer; ++i) {
    out.push_back("uptake_" + plant::util::to_string(i + 1));
  }
  return out;
}

// The whole block Jacobian, one exact column per input, by forward tangent.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_block_jacobian_forward_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                       int node) {
  const patch_type& patch = *obj_;
  tangent_block block(patch, locate(patch, static_cast<size_t>(node - 1)));
  const size_t n_in = block.in.size();
  Rcpp::NumericMatrix out(static_cast<int>(block.n_out),
                          static_cast<int>(n_in));
  for (size_t c = 0; c < n_in; ++c) {
    const std::vector<double> column = block.column(static_cast<int>(c));
    for (size_t r = 0; r < block.n_out; ++r) {
      out(static_cast<int>(r), static_cast<int>(c)) = column[r];
    }
  }
  return out;
}

// The block's value, for the non-vacuity checks: a Jacobian of a block that
// returned nothing is not evidence of anything.
// [[Rcpp::export]]
std::vector<double> ladder_block_value_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_, int node) {
  const patch_type& patch = *obj_;
  tangent_block block(patch, locate(patch, static_cast<size_t>(node - 1)));
  return block.value();
}

// The same Jacobian one row at a time, by recording the block and sweeping it
// with a unit output adjoint. This is the recording the sweep makes, so a
// disagreement with the tangent above is a disagreement about the block itself
// rather than about any transpose.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_block_jacobian_reverse_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                       int node) {
  patch_type& patch = *obj_;
  const node_address at = locate(patch, static_cast<size_t>(node - 1));

  using scalar = odelia::ode::active_scalar<double>;
  using active_strategy = plant::at_scalar<strategy_type, scalar>;
  using active_environment = plant::at_scalar<environment_type, scalar>;
  using active_individual = plant::Individual<active_strategy,
                                              active_environment>;

  const active_strategy strategy_template =
    patch.at_species(at.species).strategy_ptr()->rebind_from<scalar>();
  active_environment environment_template =
    patch.r_environment().rebind_from<scalar>();
  environment_template.light_availability.spline.set_nodes(
    patch.r_environment().light_availability.spline.knots());

  typename active_strategy::ptr sizing =
    std::make_shared<active_strategy>(strategy_template);
  active_individual probe(sizing);
  const size_t n_in = probe.block_input_size(environment_template);
  const size_t n_out = probe.block_output_size(environment_template);

  std::vector<double> in(n_in);
  patch.at_species(at.species).node_at(at.node).individual
    .block_inputs(in.begin(), patch.r_environment());

  typename scalar::tape_type tape(false);
  std::vector<std::vector<double>> rows;
  rows.reserve(n_out);
  for (size_t r = 0; r < n_out; ++r) {
    typename active_strategy::ptr strategy =
      std::make_shared<active_strategy>(strategy_template);
    active_individual individual(strategy);
    active_environment environment = environment_template;

    // A batch of one, deliberately: this reference takes the Jacobian a row at a
    // time so that what it checks is one recording swept by one unit adjoint.
    odelia::ode::adjoint_rows seed(1, n_out);
    seed[0][r] = 1.0;
    odelia::ode::adjoint_rows row;
    auto block = [&](const std::vector<scalar>& x,
                     std::vector<scalar>& y) -> void {
      individual.set_block_inputs(x.begin(), environment);
      individual.compute_rates(environment);
      individual.block_outputs(y.begin(), environment);
    };
    odelia::ode::vector_jacobian_product(tape, in, seed, block, row);
    rows.emplace_back(row[0].begin(), row[0].end());
  }
  return to_matrix(rows, n_in);
}

// The rates the reference itself computes, at the state the patch holds. A
// reference whose value disagrees with the model is a reference to a different
// function, and every Jacobian taken from it is that function's, so this is the
// first thing to compare and it needs no derivative at all.
// [[Rcpp::export]]
std::vector<double> ladder_rhs_value_forward_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  const size_t n = patch.ode_size();
  const double time = patch.ode_time();
  std::vector<double> x0(n);
  patch.ode_state(x0.begin());

  auto active = patch.rebind_from<tangent>();
  std::vector<tangent> x(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = x0[i];
  }
  active.set_ode_state(x.begin(), time);
  std::vector<tangent> dydt(n);
  active.ode_rates(dydt.begin());

  std::vector<double> out(n);
  for (size_t i = 0; i < n; ++i) {
    out[i] = odelia::util::to_passive(dydt[i]);
  }
  return out;
}

// One right-hand-side evaluation's Jacobian with respect to the ODE state, by
// forward tangent. The reference for the reduction transposes.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_rhs_state_jacobian_forward_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  const size_t n = patch.ode_size();
  const double time = patch.ode_time();
  std::vector<double> x0(n);
  patch.ode_state(x0.begin());

  std::vector<std::vector<double>> columns;
  columns.reserve(n);
  for (size_t c = 0; c < n; ++c) {
    auto active = patch.rebind_from<tangent>();
    std::vector<tangent> x(n);
    for (size_t i = 0; i < n; ++i) {
      x[i] = x0[i];
    }
    seed_direction(x[c], 1.0);
    active.set_ode_state(x.begin(), time);
    std::vector<tangent> dydt(n);
    active.ode_rates(dydt.begin());
    std::vector<double> column(n);
    for (size_t i = 0; i < n; ++i) {
      column[i] = derivative_along(dydt[i]);
    }
    columns.push_back(column);
  }
  // Columns were formed one per input; the Jacobian is rates by state.
  Rcpp::NumericMatrix out(static_cast<int>(n), static_cast<int>(n));
  for (size_t c = 0; c < n; ++c) {
    for (size_t r = 0; r < n; ++r) {
      out(static_cast<int>(r), static_cast<int>(c)) = columns[c][r];
    }
  }
  return out;
}

// The same evaluation's Jacobian with respect to the traits, species-major in
// each strategy's own order. The parameters are seeded before the state is set,
// because the quantities a state determines read the parameters and would
// otherwise be derived at the previous values.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_rhs_trait_jacobian_forward_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  const size_t n = patch.ode_size();
  const double time = patch.ode_time();
  std::vector<double> x0(n);
  patch.ode_state(x0.begin());

  size_t n_trait = 0;
  for (size_t i = 0; i < patch.size(); ++i) {
    n_trait += patch.at_species(i).strategy_ptr()->ad_parameters().size();
  }

  std::vector<std::vector<double>> columns;
  columns.reserve(n_trait);
  for (size_t c = 0; c < n_trait; ++c) {
    auto active = patch.rebind_from<tangent>();
    size_t at = 0;
    for (size_t i = 0; i < active.size(); ++i) {
      std::vector<tangent*> pars =
        active.at_species(i).strategy_ptr()->ad_parameters();
      for (size_t p = 0; p < pars.size(); ++p, ++at) {
        if (at == c) {
          seed_direction(*pars[p], 1.0);
        }
      }
    }
    std::vector<tangent> x(n);
    for (size_t i = 0; i < n; ++i) {
      x[i] = x0[i];
    }
    active.set_ode_state(x.begin(), time);
    std::vector<tangent> dydt(n);
    active.ode_rates(dydt.begin());
    std::vector<double> column(n);
    for (size_t i = 0; i < n; ++i) {
      column[i] = derivative_along(dydt[i]);
    }
    columns.push_back(column);
  }

  Rcpp::NumericMatrix out(static_cast<int>(n), static_cast<int>(n_trait));
  for (size_t c = 0; c < n_trait; ++c) {
    for (size_t r = 0; r < n; ++r) {
      out(static_cast<int>(r), static_cast<int>(c)) = columns[c][r];
    }
  }
  return out;
}

// The transpose of one right-hand-side evaluation, as the sweep takes it: the
// adjoints of the rates in, the adjoints of the state and of the traits out.
// [[Rcpp::export]]
Rcpp::List ladder_rhs_adjoint_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                   std::vector<double> lambda_dydt) {
  patch_type& patch = *obj_;
  // The flag latches for as long as the storage behind it lives, so it is cleared
  // where this evaluation starts rather than trusted to be clean.
  for (size_t i = 0; i < patch.size(); ++i) {
    const auto s = patch.at_species(i).strategy_ptr();
    *s->uptake_rows_unavailable = false;
    s->uptake_rows_reason->clear();
  }
  odelia::ode::adjoint_tape<double> tape(false);
  std::vector<double> lambda_y, trait_adjoint;
  const size_t recording =
    rhs_adjoint(patch, tape, lambda_dydt, lambda_y, trait_adjoint);
  // A recording that left its water rows off the tape returns zeros for them, and
  // a zero row and a severed one are the same number -- which is the one thing
  // this whole channel exists to prevent. So the numbers go not-a-number, the way
  // a refused metric's do, and the reason travels beside them.
  bool refused = false;
  std::string reason;
  for (size_t i = 0; i < patch.size(); ++i) {
    const auto s = patch.at_species(i).strategy_ptr();
    if (*s->uptake_rows_unavailable) {
      refused = true;
      reason = *s->uptake_rows_reason;
      break;
    }
  }
  if (refused) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    lambda_y.assign(lambda_y.size(), nan);
    trait_adjoint.assign(trait_adjoint.size(), nan);
  }
  // The knot halves are gone with the reduction transposes that produced them:
  // the field's rows are an intermediate of the stage recording now, so they are
  // in the state and trait rows rather than beside them.
  return Rcpp::List::create(Rcpp::_["state"] = lambda_y,
                            Rcpp::_["trait"] = trait_adjoint,
                            Rcpp::_["refused"] = refused,
                            Rcpp::_["reason"] = reason,
                            Rcpp::_["block_recording_size"] =
                              static_cast<double>(recording));
}

// The trajectory reference: one exact directional derivative of the census by a
// forward tangent of the same run, stepped at the sizes it recorded. `direction`
// carries one weight per trait column, so a unit vector gives one Jacobian column
// and a mixed direction gives the contraction the trajectory rungs compare.
// [[Rcpp::export]]
Rcpp::List ladder_trajectory_tangent_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                          std::vector<double> direction) {
  std::vector<double> value;
  const std::vector<double> tangent =
    obj_->census_trait_tangent(direction, value);
  return Rcpp::List::create(Rcpp::_["value"] = value,
                            Rcpp::_["tangent"] = tangent);
}

// The inflow boundary's density at an active scalar, with its row in one
// registered parameter, read where the condition forms it. `index` is a position
// in ad_parameters(), counting from one.
//
// The boundary node is not ODE state, so nothing downstream reports what the
// differentiated path made of this row; a caller comparing it against a rebuilt
// difference of the same quantity is comparing the row against the condition.
// [[Rcpp::export]]
Rcpp::List ladder_boundary_density_tangent_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                int index) {
    const patch_type& source = *obj_;
  auto active = source.template rebind_from<tangent>();
  std::vector<tangent*> pars =
    active.at_species(0).strategy_ptr()->ad_parameters();
  if (index < 1 || static_cast<size_t>(index) > pars.size()) {
    plant::util::stop("ladder_boundary_density_tangent: index out of range");
  }
  seed_direction(*pars[static_cast<size_t>(index) - 1], 1.0);

  // Seeded before the state is applied, because an active patch derives the aux a state
  // determines as it is built, which is before this seed exists. Leaf area at
  // fixed height would then carry no row, and neither would the field the
  // reduction builds from it -- so the two constants that set leaf area come back
  // short by their whole field channel while every other parameter is exact.
  const size_t n = source.ode_size();
  std::vector<double> x0(n);
  source.ode_state(x0.begin());
  std::vector<tangent> x(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = x0[i];
  }
  active.set_ode_state(x.begin(), source.ode_time());

  // A rate evaluation owns the field build and the boundary condition, so this
  // leaves the node where a step would have left it.
  std::vector<tangent> dydt(active.ode_size());
  active.r_compute_environment();
  active.ode_rates(dydt.begin());

  const auto& node = active.at_species(0).r_new_node();
  const tangent ell = node.get_log_density();
  const tangent h = node.height();
  // The birth-size carbon the establishment probability divides by, read where
  // the condition reads it.
  const tangent carbon = node.individual.aux("net_mass_production_dt");
  const tangent area = node.individual.aux("competition_effect");
  return Rcpp::List::create(
    Rcpp::_["log_density"] = odelia::util::to_passive(ell),
    Rcpp::_["dlog_density"] = derivative_along(ell),
    Rcpp::_["height"] = odelia::util::to_passive(h),
    Rcpp::_["dheight"] = derivative_along(h),
    Rcpp::_["carbon"] = odelia::util::to_passive(carbon),
    Rcpp::_["dcarbon"] = derivative_along(carbon),
    Rcpp::_["area_leaf"] = odelia::util::to_passive(area),
    Rcpp::_["darea_leaf"] = derivative_along(area));
}

// The seed's height and leaf area at an active scalar, with their row in one
// registered parameter, read where the recording forms it rather than inferred from a
// census. `index` is a position in ad_parameters(), counting from one.
//
// The referee is a difference of the seed height over a rebuilt strategy: the recording
// declares this row, so a caller comparing the two is comparing a declaration
// against the condition it claims to solve.
// [[Rcpp::export]]
Rcpp::List ladder_seed_geometry_tangent_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                             int index) {
    const patch_type& patch = *obj_;
  const strategy_type& source = *patch.at_species(0).strategy_ptr();
  auto active = source.template rebind_from<tangent>();
  std::vector<tangent*> pars = active.ad_parameters();
  if (index < 1 || static_cast<size_t>(index) > pars.size()) {
    plant::util::stop("ladder_seed_geometry_tangent: index out of range");
  }
  seed_direction(*pars[static_cast<size_t>(index) - 1], 1.0);
  const auto seed = active.seed_geometry();
  return Rcpp::List::create(
    Rcpp::_["height"] = odelia::util::to_passive(seed.height),
    Rcpp::_["dheight"] = derivative_along(seed.height),
    Rcpp::_["area_leaf"] = odelia::util::to_passive(seed.area_leaf),
    Rcpp::_["darea_leaf"] = derivative_along(seed.area_leaf));
}

// One exact directional derivative of the census with respect to the first
// recorded state. `direction` carries one weight per component of that state, so a
// coordinate direction gives the census's sensitivity to the value one cohort
// starts at, with no trait and no derived quantity on the path.
// [[Rcpp::export]]
Rcpp::List ladder_census_initial_state_tangent_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                    std::vector<double> direction,
                                                    int segment) {
  std::vector<double> value;
  const std::vector<double> tangent =
    obj_->census_initial_state_tangent(
      direction, value, static_cast<size_t>(segment));
  return Rcpp::List::create(Rcpp::_["value"] = value,
                            Rcpp::_["tangent"] = tangent);
}

// The state a segment's first step ran from, which is what `direction` and
// `state0` are indexed against. A widened state is not what any record holds, so
// a caller cannot read it off the trajectory.
// [[Rcpp::export]]
std::vector<double> ladder_segment_base_state_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                   int segment) {
  return obj_->segment_base_state(static_cast<size_t>(segment));
}

// The census a plain-double replay of the recorded steps reaches from `state0`.
// Differencing this is what referees the tangent above: the same steps, the same
// introductions, one perturbed state.
// [[Rcpp::export]]
std::vector<double> ladder_census_initial_state_replay_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                            std::vector<double> state0,
                                                            int segment) {
  return obj_->census_initial_state_replay(
    state0, static_cast<size_t>(segment));
}

// How many times the inflow boundary's own term entered the trait adjoint over one
// sweep, beside the number of metrics that sweep produced.
//
// The boundary node stands at the seed's height for the whole run and its condition
// is re-evaluated at every stage of every step, so this row is multiplied by a
// count no other check reads. A row correct per evaluation and wrong in what it is
// multiplied by is a different failure from a wrong row, and neither a forward
// tangent nor a sweep can see it, because both apply the same multiplier.
//
// The condition sits inside the stage recording, so the count is the solver's
// count of stage transposes: every one records it and every sweep of one carries
// its row.
// [[Rcpp::export]]
Rcpp::List ladder_boundary_evaluations_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  obj_->clear_boundary_condition_evaluations();
  const plant::census_gradient result =
    obj_->census_trait_gradient();
  // And how many operating points the sweep placed rather than searched for. A
  // record that engages and one that quietly does not produce the same numbers, so
  // this is the only thing that tells them apart.
  double placed = 0.0;
  for (size_t i = 0; i < obj_->r_patch().size(); ++i) {
    placed += static_cast<double>(
      obj_->r_patch().at_species(i).strategy_ptr()->leaf_placements());
  }
  return Rcpp::List::create(
    Rcpp::_["evaluations"] =
      static_cast<double>(obj_->boundary_condition_evaluations()),
    Rcpp::_["placements"] = placed,
    Rcpp::_["metrics"] = static_cast<int>(result.gradient.size()));
}

// The trait columns' names, species-major and each carrying its species index,
// so a gradient's columns can be told apart when two species carry the same
// parameter.
// [[Rcpp::export]]
std::vector<std::string> ladder_trait_names_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->trait_adjoint_names();
}

// How many nodes the block loop visits, which is what a node index runs over.
// [[Rcpp::export]]
int ladder_node_count_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return static_cast<int>(obj_->node_count());
}

// The field's knot data, which is what the recorded step reads. The reduction
// that builds it has no transpose of its own any more -- it is an intermediate
// of the stage recording -- so what this serves is the forward check that a
// permutation of the nodes leaves the knots alone.
// [[Rcpp::export]]
Rcpp::List ladder_field_knots_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  const environment_type& environment = patch.r_environment();
  const std::vector<double>& value = environment.light_availability.knot_values();
  const std::vector<double>& slope = environment.light_availability.knot_slopes();
  return Rcpp::List::create(Rcpp::_["value"] = value,
                            Rcpp::_["slope"] = slope,
                            Rcpp::_["height"] =
                              environment.light_availability.spline.knots());
}


// The block's outputs differenced in one input, in plain double, with the
// environment held. A difference of the RECORDED step cannot see an error in a
// supplied row -- the recording makes the value independent of a recorded input -- so
// this differences the double path instead, which runs the leaf's own solve.
// Columns are the block's inputs in block_inputs order, rows its outputs.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_block_difference_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                 int node, double rel) {
  patch_type& patch = *obj_;
  const node_address at = locate(patch, static_cast<size_t>(node - 1));

  const strategy_type& source = *patch.at_species(at.species).strategy_ptr();
  environment_type environment = patch.r_environment();
  strategy_type::ptr strategy = std::make_shared<strategy_type>(source);
  plant::Individual<strategy_type, environment_type> individual(strategy);

  const size_t n_in = individual.block_input_size(environment);
  const size_t n_out = individual.block_output_size(environment);
  std::vector<double> in(n_in);
  patch.at_species(at.species).node_at(at.node).individual
    .block_inputs(in.begin(), patch.r_environment());

  auto evaluate = [&](const std::vector<double>& x, std::vector<double>& y) {
    environment = patch.r_environment();
    individual.set_block_inputs(x.begin(), environment);
    individual.compute_rates(environment);
    y.resize(n_out);
    individual.block_outputs(y.begin(), environment);
  };

  Rcpp::NumericMatrix out(static_cast<int>(n_out), static_cast<int>(n_in));
  std::vector<double> up(n_out), dn(n_out);
  for (size_t c = 0; c < n_in; ++c) {
    const double h = std::max(std::abs(in[c]) * rel, rel);
    std::vector<double> x = in;
    x[c] = in[c] + h;
    evaluate(x, up);
    x[c] = in[c] - h;
    evaluate(x, dn);
    for (size_t r = 0; r < n_out; ++r) {
      out(static_cast<int>(r), static_cast<int>(c)) = (up[r] - dn[r]) / (2.0 * h);
    }
  }
  return out;
}

// The block's outputs differenced along ONE direction of its inputs, in plain
// double, with the environment held.
//
// Summing the per-column differences above would not do. Water moves on
// differences of potential and tissue fails on absolutes, so along a direction
// that moves every soil potential together the answer is a small residue of the
// columns that make it -- and a sum of columns carries the columns' own error,
// which is the size of the terms rather than the size of what is left of them.
// Perturbing along the direction puts that cancellation inside the model, where
// the same near-symmetry shrinks the truncation too.
//
// `direction` is one weight per block input, in block_inputs order; the step is
// `rel` times the size of the inputs the direction actually reaches.
// [[Rcpp::export]]
std::vector<double> ladder_block_direction_difference_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                           int node,
                                                           std::vector<double> direction,
                                                           double rel) {
  patch_type& patch = *obj_;
  const node_address at = locate(patch, static_cast<size_t>(node - 1));

  const strategy_type& source = *patch.at_species(at.species).strategy_ptr();
  environment_type environment = patch.r_environment();
  strategy_type::ptr strategy = std::make_shared<strategy_type>(source);
  plant::Individual<strategy_type, environment_type> individual(strategy);

  const size_t n_in = individual.block_input_size(environment);
  const size_t n_out = individual.block_output_size(environment);
  plant::util::check_length(direction.size(), n_in);
  std::vector<double> in(n_in);
  patch.at_species(at.species).node_at(at.node).individual
    .block_inputs(in.begin(), patch.r_environment());

  // The step is scaled by the inputs the direction reaches, so a direction over
  // one family is not stepped by the size of another.
  double reach = 0.0;
  for (size_t c = 0; c < n_in; ++c) {
    if (direction[c] != 0.0) {
      reach = std::max(reach, std::abs(in[c]));
    }
  }
  if (reach == 0.0) {
    plant::util::stop("a direction reaching no input of any size");
  }
  const double h = reach * rel;

  auto evaluate = [&](const std::vector<double>& x, std::vector<double>& y) {
    environment = patch.r_environment();
    individual.set_block_inputs(x.begin(), environment);
    individual.compute_rates(environment);
    y.resize(n_out);
    individual.block_outputs(y.begin(), environment);
  };

  std::vector<double> up(n_out), dn(n_out), x(n_in);
  for (size_t c = 0; c < n_in; ++c) {
    x[c] = in[c] + h * direction[c];
  }
  evaluate(x, up);
  for (size_t c = 0; c < n_in; ++c) {
    x[c] = in[c] - h * direction[c];
  }
  evaluate(x, dn);

  std::vector<double> out(n_out);
  for (size_t r = 0; r < n_out; ++r) {
    out[r] = (up[r] - dn[r]) / (2.0 * h);
  }
  return out;
}

// The introduction map's whole Jacobian, both ways, at one widening.
//
// This is rung 5's own unit and the one object it had no reference for. The map is
// the pre-introduction state and the traits in, the whole widened state out; the
// forward side seeds one tangent per input column, and the reverse side is the
// transpose the sweep actually runs, seeded with one unit output adjoint per row.
// Both go through introduce_over, so the reference traverses the forward function
// and not the transpose.
//
// Forming it entirely is what localises a disagreement to a cell. A contraction
// returns one number, hides an error behind a small seed component, and when it
// fails localises to nothing.
// [[Rcpp::export]]
Rcpp::List ladder_introduction_jacobian_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                             std::vector<int> species_index,
                                             std::vector<double> state_before,
                                             double time_before) {
  patch_type& patch = *obj_;
  std::vector<size_t> which;
  which.reserve(species_index.size());
  for (const int i : species_index) {
    if (i < 1 || static_cast<size_t>(i) > patch.size()) {
      plant::util::stop("a species index is counted from one");
    }
    which.push_back(static_cast<size_t>(i - 1));
  }

  const size_t n_state = patch.ode_size();
  const size_t n_trait = patch.trait_adjoint_size();
  const std::vector<std::vector<double>> forward =
    patch.introduction_jacobian(which, state_before, time_before);
  const size_t n_out = forward.size();

  // One seed per output row, swept over a single recording of the map, which is
  // the shape the trajectory sweep takes. The trait accumulator adds by design,
  // so it is cleared before the call and again after it.
  const odelia::ode::adjoint_rows seeds = odelia::ode::adjoint_rows::all_rows(n_out);
  // The transpose the sweep takes at a widening, taken here on its own: the same
  // tape and the same product the solver's walk runs, with the trajectory and the
  // segmentation left out so a disagreement is the map's.
  using ad_scalar = odelia::ode::active_scalar<double>;
  auto widen = [&](auto& active_patch,
                   std::vector<ad_scalar>::const_iterator x,
                   std::vector<ad_scalar>& y) -> void {
    active_patch.apply_insertion(which, time_before, x, y);
  };
  odelia::ode::adjoint_rows lambda_before;
  odelia::ode::adjoint_rows trait_adjoint(n_out, n_trait);
  ad_scalar::tape_type tape(false);
  odelia::ode::active_system<std::decay_t<decltype(patch)>> active{patch, tape};
  odelia::ode::state_and_parameter_adjoints(active, state_before, seeds, widen,
                                            lambda_before, trait_adjoint);

  Rcpp::NumericMatrix rev(static_cast<int>(n_out),
                          static_cast<int>(n_state + n_trait));
  for (size_t r = 0; r < n_out; ++r) {
    for (size_t c = 0; c < n_state; ++c) {
      rev(static_cast<int>(r), static_cast<int>(c)) = lambda_before[r][c];
    }
    for (size_t p = 0; p < n_trait; ++p) {
      rev(static_cast<int>(r), static_cast<int>(n_state + p)) =
        trait_adjoint[r][p];
    }
  }

  Rcpp::NumericMatrix fwd(static_cast<int>(n_out),
                          static_cast<int>(n_state + n_trait));
  for (size_t r = 0; r < n_out; ++r) {
    for (size_t c = 0; c < n_state + n_trait; ++c) {
      fwd(static_cast<int>(r), static_cast<int>(c)) = forward[r][c];
    }
  }
  return Rcpp::List::create(Rcpp::_["forward"] = fwd, Rcpp::_["reverse"] = rev,
                            Rcpp::_["traits"] = patch.trait_adjoint_names());
}

// The whole right-hand side differenced in each trait, in plain double, by
// perturbing the prepared strategy exactly where the tangent seeds it. A
// difference that rebuilds from Parameters re-runs prepare_strategy and so
// carries the birth-size channel the differentiated path imposes to zero; this
// one does not, so it referees the trait rows the sweep actually computes.
// [[Rcpp::export]]
Rcpp::NumericMatrix ladder_rhs_trait_difference_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                                     double rel) {
  patch_type& patch = *obj_;
  const size_t n = patch.ode_size();
  const double time = patch.ode_time();
  std::vector<double> x0(n);
  patch.ode_state(x0.begin());

  std::vector<double*> pars;
  for (size_t i = 0; i < patch.size(); ++i) {
    for (double* p : patch.at_species(i).strategy_ptr()->ad_parameters()) {
      pars.push_back(p);
    }
  }

  auto rates_at = [&](std::vector<double>& y) {
    patch.set_ode_state(x0.begin(), time);
    y.resize(n);
    patch.ode_rates(y.begin());
  };

  Rcpp::NumericMatrix out(static_cast<int>(n), static_cast<int>(pars.size()));
  std::vector<double> up(n), dn(n);
  for (size_t c = 0; c < pars.size(); ++c) {
    const double base = *pars[c];
    const double h = std::max(std::abs(base) * rel, rel);
    *pars[c] = base + h;
    rates_at(up);
    *pars[c] = base - h;
    rates_at(dn);
    *pars[c] = base;
    for (size_t r = 0; r < n; ++r) {
      out(static_cast<int>(r), static_cast<int>(c)) = (up[r] - dn[r]) / (2.0 * h);
    }
  }
  rates_at(up);
  return out;
}


// A rebind and an assignment have to leave a patch holding the same thing: they
// are one map reached two ways, and every comment in the tree says so. Nothing
// checked it, and they had already drifted.
//
// A rebind hands back a fresh object, so a member it does not write is
// default-constructed; an assignment writes into one that exists, so a member it
// does not write keeps what it had. Every member is a place they can differ.
//
// What is compared is not what the two arrive holding. Neither is required to
// arrive with a field: the light spline is derived from the state, and both
// operations leave building it to the caller. What they are required to do is
// reach the same place from the same state, so both are put through the call
// every caller makes -- and then the rates are compared, because the rates are
// what a stage recording records.
// [[Rcpp::export]]
Rcpp::List ladder_rebind_matches_assign_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;

  std::vector<double> y(patch.ode_size());
  patch.ode_state(y.begin());
  const double t = patch.time();

  adjoint_patch rebound = patch.template rebind_from<adjoint>();

  // Assigned into a patch that has already been used, because that is the only
  // patch the sweep ever assigns into: one holding the last recording's state,
  // its field, and scalars carrying that recording's tape slots. Assigning into a
  // fresh one tests the easy half and passes while the real path is wrong.
  adjoint_patch assigned = patch.template rebind_from<adjoint>();
  {
    std::vector<adjoint> used(y.size());
    for (size_t i = 0; i < y.size(); ++i) { used[i] = adjoint(y[i] * 1.01); }
    assigned.set_ode_state(used.begin(), t + 1.0);
    std::vector<adjoint> scratch(assigned.ode_size());
    assigned.ode_rates(scratch.begin());
  }
  assigned.assign_from(patch);

  std::vector<adjoint> x(y.size());
  for (size_t i = 0; i < y.size(); ++i) { x[i] = adjoint(y[i]); }
  rebound.set_ode_state(x.begin(), t);
  assigned.set_ode_state(x.begin(), t);

  auto values = [](const std::vector<adjoint>& v) {
    std::vector<double> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
      out[i] = odelia::util::to_passive(v[i]);
    }
    return out;
  };
  auto state_of = [&](adjoint_patch& p) {
    std::vector<adjoint> raw(p.ode_size());
    p.ode_state(raw.begin());
    return values(raw);
  };
  auto rates_of = [&](adjoint_patch& p) {
    std::vector<adjoint> raw(p.ode_size());
    p.ode_rates(raw.begin());
    return values(raw);
  };
  auto field_of = [&](adjoint_patch& p) {
    const auto env = p.r_environment();
    std::vector<adjoint> raw(env.n_cohort_reads());
    // Held against the width rather than discarded: the buffer starts at zero, so
    // a fill that writes nothing reads the same as one that writes zeros, and the
    // comparison this feeds would pass on two buffers neither side had filled.
    const auto end = env.cohort_reads(raw.begin());
    plant::util::check_length(
      static_cast<size_t>(std::distance(raw.begin(), end)), raw.size());
    return values(raw);
  };
  auto parameters_of = [](adjoint_patch& p) {
    std::vector<double> out;
    for (const adjoint* q : p.ad_parameters()) {
      out.push_back(odelia::util::to_passive(*q));
    }
    return out;
  };
  auto worst = [](const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
      return std::numeric_limits<double>::infinity();
    }
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      const double d = std::abs(a[i] - b[i]);
      if (d > m) { m = d; }
    }
    return m;
  };

  return Rcpp::List::create(
    Rcpp::_["ode_size"] = Rcpp::IntegerVector::create((int) rebound.ode_size(),
                                                      (int) assigned.ode_size()),
    Rcpp::_["aux_size"] = Rcpp::IntegerVector::create((int) rebound.aux_size(),
                                                      (int) assigned.aux_size()),
    Rcpp::_["n_cohort_reads"] =
      Rcpp::IntegerVector::create((int) rebound.r_environment().n_cohort_reads(),
                                  (int) assigned.r_environment().n_cohort_reads()),
    Rcpp::_["state_gap"] = worst(state_of(rebound), state_of(assigned)),
    Rcpp::_["rate_gap"] = worst(rates_of(rebound), rates_of(assigned)),
    Rcpp::_["field_gap"] = worst(field_of(rebound), field_of(assigned)),
    Rcpp::_["parameter_gap"] = worst(parameters_of(rebound),
                                     parameters_of(assigned)));
}

// Defined in census_gradient.cpp beside the whole-sweep entry point, so a
// split sweep and a whole one reach R in one shape.
Rcpp::List census_gradient_to_r(const plant::census_gradient& g);

// The census's own reading of the traits at the state held. No sweep produces it,
// so it is the one route to the census a transpose check cannot touch, and it is
// easy to omit because it is a one-line calculation at the final state.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_direct_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_state_and_trait_rows().trait.to_rows();
}

// The same quantity differenced in plain double, which is what referees it.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_difference_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                             double rel) {
  return obj_->census_trait_difference(rel);
}

// The same gradient with the sweep stopped and resumed at each given recorded
// step, counted from one. Composition over steps is associative, so this must
// agree with census_trait_gradient_tf24 bit for bit.
// [[Rcpp::export]]
Rcpp::List
census_trait_gradient_split_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                 std::vector<int> splits) {
  std::vector<size_t> at;
  at.reserve(splits.size());
  for (const int s : splits) {
    if (s < 1) {
      plant::util::stop("a split point is a recorded step counted from one");
    }
    at.push_back(static_cast<size_t>(s - 1));
  }
  return census_gradient_to_r(obj_->census_trait_gradient(at));
}

// How many backward ranges the last gradient swept. A requested split that fell
// on a segment boundary rather than inside one changes nothing, so this is what
// says a split actually cut.
// [[Rcpp::export]]
double census_adjoint_segments_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return static_cast<double>(obj_->adjoint_segments);
}

// The adjoint the last gradient's walk ended holding: d(census)/d(the first
// recorded state), one row per metric swept and one column per entry of that
// state. Reaching it means carrying lambda over the range below the first
// widening, so a run with steps there is what separates a walk that ran that
// range from one that started above it.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_adjoint_at_first_state_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->adjoint_at_first_state;
}

// The forward run's classification tally: one row per species, one column per
// operating-point kind in census_operating_point_names_tf24() order.
//
// The leaf classifies by the branch taken and the next plant overwrites it, so a
// run's incidence is not recoverable afterwards from anything but this. It is
// what says whether a regime the gradient refuses is rare or is most of the run.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_operating_point_counts_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const std::vector<std::vector<size_t>> counts = obj_->operating_point_counts();
  std::vector<std::vector<double>> ret;
  ret.reserve(counts.size());
  for (const std::vector<size_t>& row : counts) {
    ret.push_back(std::vector<double>(row.begin(), row.end()));
  }
  return ret;
}

// [[Rcpp::export]]
void census_clear_operating_point_counts_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  obj_->clear_operating_point_counts();
}

// The kinds, in the order the counts are reported.
// [[Rcpp::export]]
std::vector<std::string> census_operating_point_names_tf24() {
  std::vector<std::string> ret;
  ret.reserve(plant::Leaf::operating_point_kind_count);
  for (size_t k = 0; k < plant::Leaf::operating_point_kind_count; ++k) {
    ret.push_back(plant::Leaf::operating_point_kind_name(
        static_cast<plant::Leaf::OperatingPointKind>(k)));
  }
  return ret;
}

// How often each counted clamp fired on the forward run, one row per species.
// A clamp masking a smooth function severs a gradient row, and a severed row and
// a true zero are the same number -- so the honest treatment is to count it.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_clamp_counts_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const std::vector<std::vector<size_t>> counts = obj_->clamp_counts();
  std::vector<std::vector<double>> ret;
  ret.reserve(counts.size());
  for (const std::vector<size_t>& row : counts) {
    ret.push_back(std::vector<double>(row.begin(), row.end()));
  }
  return ret;
}

// The same sites counted where the sweep runs. A clamp only severs a row on the
// differentiated path, so this is the tally that says whether a gradient carries
// a severance -- the forward one says only that the guard is reachable.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_clamp_counts_differentiated_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const std::vector<std::vector<size_t>> counts = obj_->clamp_counts_differentiated();
  std::vector<std::vector<double>> ret;
  ret.reserve(counts.size());
  for (const std::vector<size_t>& row : counts) {
    ret.push_back(std::vector<double>(row.begin(), row.end()));
  }
  return ret;
}

// The clamp sites, in the order the counts are reported. Read from the enum so a
// site cannot be counted under its neighbour's name.
// [[Rcpp::export]]
std::vector<std::string> census_clamp_names_tf24() {
  std::vector<std::string> ret;
  ret.reserve(plant::CLAMP_SITE_COUNT);
  for (int i = 0; i < plant::CLAMP_SITE_COUNT; ++i) {
    ret.push_back(plant::clamp_site_name(i));
  }
  return ret;
}

// The smallest profit curvature the differentiated path met, one per species, or
// -1 where it met none. The guard on it refuses a row rather than returning
// amplification, and a guard that held reports the same green as a guard nothing
// reached -- so the distance to the floor is reported rather than assumed.
// [[Rcpp::export]]
std::vector<double>
census_curvature_margin_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->curvature_margins();
}
