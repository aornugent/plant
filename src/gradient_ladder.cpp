#include <plant.h>
#include <XAD/XAD.hpp>
#include <chrono>

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
// referees is the graft. Taking it forward and backward checks that both modes
// traverse the same recording, and it cannot check the number a grafted row
// carries: both modes read that number from the same place, so they agree
// whether it is right or wrong. The plain-double difference below is the only
// instrument here that can see a wrong supplied row, because it runs the leaf's
// own solve instead of reading the row.

namespace {

using strategy_type = plant::TF24_Strategy<double>;
using environment_type = plant::TF24_Environment<double>;
using patch_type = plant::Patch<strategy_type, environment_type>;
using patch_handle = plant::RcppR6::RcppR6<patch_type>;

using tangent = xad::fwd<double>::active_type;
using tangent_strategy = strategy_type::rebind<tangent>;
using tangent_environment = environment_type::rebind<tangent>;
using tangent_individual = plant::Individual<tangent_strategy,
                                             tangent_environment>;

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
      xad::derivative(x[static_cast<size_t>(seed)]) = 1.0;
    }

    individual.set_block_inputs(x.begin(), environment);
    individual.compute_rates(environment);
    std::vector<tangent> y(n_out);
    individual.block_outputs(y.begin(), environment);

    std::vector<double> out(n_out);
    for (size_t j = 0; j < n_out; ++j) {
      out[j] = xad::derivative(y[j]);
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
      out[j] = xad::value(y[j]);
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
  using active_strategy = strategy_type::rebind<scalar>;
  using active_environment = environment_type::rebind<scalar>;
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

    std::vector<double> seed(n_out, 0.0), row;
    seed[r] = 1.0;
    auto block = [&](const std::vector<scalar>& x,
                     std::vector<scalar>& y) -> void {
      individual.set_block_inputs(x.begin(), environment);
      individual.compute_rates(environment);
      individual.block_outputs(y.begin(), environment);
    };
    odelia::ode::vector_jacobian_product(tape, in, seed, block, row);
    rows.push_back(row);
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
    out[i] = xad::value(dydt[i]);
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
    xad::derivative(x[c]) = 1.0;
    active.set_ode_state(x.begin(), time);
    std::vector<tangent> dydt(n);
    active.ode_rates(dydt.begin());
    std::vector<double> column(n);
    for (size_t i = 0; i < n; ++i) {
      column[i] = xad::derivative(dydt[i]);
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
          xad::derivative(*pars[p]) = 1.0;
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
      column[i] = xad::derivative(dydt[i]);
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
  plant::util::check_length(lambda_dydt.size(), patch.ode_size());
  patch.clear_trait_adjoint();
  std::vector<double> lambda_y(patch.ode_size(), 0.0);
  patch.ode_rates_adjoint(lambda_dydt.begin(), lambda_y.begin());
  // The knot halves are gone with the reduction transposes that produced them:
  // the field's rows are an intermediate of the stage recording now, so they are
  // in the state and trait rows rather than beside them.
  return Rcpp::List::create(Rcpp::_["state"] = lambda_y,
                            Rcpp::_["trait"] = patch.trait_adjoint[0],
                            Rcpp::_["block_recording_size"] =
                              static_cast<double>(patch.block_recording_size),
                            Rcpp::_["block_sweeps"] =
                              static_cast<double>(patch.block_sweeps));
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
    obj_->census_trait_tangent<plant::tf24_census>(direction, value);
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
  using tangent = xad::fwd<double>::active_type;
  const patch_type& source = *obj_;
  auto active = source.template rebind_from<tangent>();
  std::vector<tangent*> pars =
    active.at_species(0).strategy_ptr()->ad_parameters();
  if (index < 1 || static_cast<size_t>(index) > pars.size()) {
    plant::util::stop("ladder_boundary_density_tangent: index out of range");
  }
  xad::derivative(*pars[static_cast<size_t>(index) - 1]) = 1.0;

  // Seeded before the state is applied, because a twin derives the aux a state
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
    Rcpp::_["log_density"] = xad::value(ell),
    Rcpp::_["dlog_density"] = xad::derivative(ell),
    Rcpp::_["height"] = xad::value(h),
    Rcpp::_["dheight"] = xad::derivative(h),
    Rcpp::_["carbon"] = xad::value(carbon),
    Rcpp::_["dcarbon"] = xad::derivative(carbon),
    Rcpp::_["area_leaf"] = xad::value(area),
    Rcpp::_["darea_leaf"] = xad::derivative(area));
}

// The seed's height and leaf area at an active scalar, with their row in one
// registered parameter, read where the graft forms it rather than inferred from a
// census. `index` is a position in ad_parameters(), counting from one.
//
// The referee is a difference of the seed height over a rebuilt strategy: the graft
// declares this row, so a caller comparing the two is comparing a declaration
// against the condition it claims to solve.
// [[Rcpp::export]]
Rcpp::List ladder_seed_geometry_tangent_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                             int index) {
  using tangent = xad::fwd<double>::active_type;
  const patch_type& patch = *obj_;
  const strategy_type& source = *patch.at_species(0).strategy_ptr();
  auto active = source.template rebind_from<tangent>();
  std::vector<tangent*> pars = active.ad_parameters();
  if (index < 1 || static_cast<size_t>(index) > pars.size()) {
    plant::util::stop("ladder_seed_geometry_tangent: index out of range");
  }
  xad::derivative(*pars[static_cast<size_t>(index) - 1]) = 1.0;
  const auto seed = active.seed_geometry();
  return Rcpp::List::create(
    Rcpp::_["height"] = xad::value(seed.height),
    Rcpp::_["dheight"] = xad::derivative(seed.height),
    Rcpp::_["area_leaf"] = xad::value(seed.area_leaf),
    Rcpp::_["darea_leaf"] = xad::derivative(seed.area_leaf));
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
    obj_->census_initial_state_tangent<plant::tf24_census>(
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
  return obj_->census_initial_state_replay<plant::tf24_census>(
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
// `asked` counts every call and `carried` the calls that recorded something; they
// differ exactly on the sweeps where no boundary adjoint was seeded.
// [[Rcpp::export]]
Rcpp::List ladder_boundary_evaluations_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  obj_->clear_boundary_condition_evaluations();
  const std::vector<std::vector<double>> gradient =
    obj_->census_trait_gradient<plant::tf24_census>();
  const std::vector<size_t> counts = obj_->boundary_condition_evaluations();
  return Rcpp::List::create(
    Rcpp::_["asked"] = static_cast<double>(counts[0]),
    Rcpp::_["carried"] = static_cast<double>(counts[1]),
    Rcpp::_["metrics"] = static_cast<int>(gradient.size()));
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
// supplied row -- the graft makes the value independent of a grafted input -- so
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
  std::vector<std::vector<double>> seeds(n_out, std::vector<double>(n_out, 0.0));
  for (size_t r = 0; r < n_out; ++r) {
    seeds[r][r] = 1.0;
  }
  std::vector<std::vector<double>> lambda_before;
  patch.clear_trait_adjoint(n_out);
  patch.introduction_adjoint(which, state_before, time_before, seeds,
                             lambda_before);

  Rcpp::NumericMatrix rev(static_cast<int>(n_out),
                          static_cast<int>(n_state + n_trait));
  for (size_t r = 0; r < n_out; ++r) {
    for (size_t c = 0; c < n_state; ++c) {
      rev(static_cast<int>(r), static_cast<int>(c)) = lambda_before[r][c];
    }
    for (size_t p = 0; p < n_trait; ++p) {
      rev(static_cast<int>(r), static_cast<int>(n_state + p)) =
        patch.trait_adjoint[r][p];
    }
  }
  patch.clear_trait_adjoint();

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

// Where one whole right-hand-side adjoint spends its time.
//
// The stage is one recording swept once per metric, so there are no components
// left to attribute between -- what used to be eight entries here were the
// hand-written transposes it replaced. What is still separable is the forward
// field build against the recording that repeats it at an active scalar, which
// is the ratio that says what differentiating the stage costs over running it.
// [[Rcpp::export]]
Rcpp::List ladder_rhs_adjoint_timing_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                          std::vector<double> lambda_dydt,
                                          int reps) {
  using clock = std::chrono::steady_clock;
  patch_type& patch = *obj_;
  plant::util::check_length(lambda_dydt.size(), patch.ode_size());

  const size_t n = patch.ode_size();
  double t_env = 0, t_total = 0;
  double slots = 0, sweeps = 0;

  for (int rep = 0; rep < reps; ++rep) {
    std::vector<double> lambda_state(n, 0.0);

    auto t0 = clock::now();
    patch.r_compute_environment();
    auto t1 = clock::now();
    t_env += std::chrono::duration<double>(t1 - t0).count();

    t0 = clock::now();
    patch.ode_rates_adjoint(lambda_dydt.begin(), lambda_state.begin());
    t1 = clock::now();
    t_total += std::chrono::duration<double>(t1 - t0).count();

    slots = static_cast<double>(patch.block_recording_size);
    sweeps = static_cast<double>(patch.block_sweeps);
  }

  const double r = static_cast<double>(reps);
  return Rcpp::List::create(
    Rcpp::_["compute_environment"] = t_env / r,
    Rcpp::_["ode_rates_adjoint"] = t_total / r,
    Rcpp::_["total"] = t_total / r,
    Rcpp::_["block_recording_size"] = slots,
    Rcpp::_["block_sweeps"] = sweeps);
}

// What one recorded block pays before it records anything.
//
// cohort_block_adjoint builds a fresh active strategy and a fresh active
// environment per cohort per stage, then registers 185 inputs and clears the
// tape. None of that scales with the recorded operation count, so a block cost
// that does not move with the quadrature rule is this rather than the tape. The
// two are timed apart here because the fix differs: one is amortisable across a
// stage, the other is not.
// [[Rcpp::export]]
Rcpp::List ladder_block_copy_cost_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                                       int node, int reps) {
  using clock = std::chrono::steady_clock;
  using ad_scalar = odelia::ode::active_scalar<double>;
  using ad_strategy = strategy_type::rebind<ad_scalar>;
  using ad_environment = environment_type::rebind<ad_scalar>;

  patch_type& patch = *obj_;
  const node_address at = locate(patch, static_cast<size_t>(node - 1));

  const ad_strategy strategy_template =
    patch.at_species(at.species).strategy_ptr()->template rebind_from<ad_scalar>();
  const ad_environment environment_template =
    patch.r_environment().template rebind_from<ad_scalar>();

  double sink = 0.0;

  // The strategy copy alone.
  auto t0 = clock::now();
  for (int i = 0; i < reps; ++i) {
    auto s = std::make_shared<ad_strategy>(strategy_template);
    sink += odelia::util::to_passive(s->pars.lma);
  }
  const double t_strategy =
    std::chrono::duration<double>(clock::now() - t0).count() / reps;

  // The environment copy alone.
  t0 = clock::now();
  for (int i = 0; i < reps; ++i) {
    ad_environment e = environment_template;
    sink += static_cast<double>(e.light_availability.spline.knots().size());
  }
  const double t_environment =
    std::chrono::duration<double>(clock::now() - t0).count() / reps;

  // A bare tape cycle over the same input count, recording nothing.
  const size_t n_in = 185;
  std::vector<double> in(n_in, 1.0), out_adjoint(12, 1.0), in_adjoint;
  typename ad_scalar::tape_type tape(false);
  auto nothing = [&](const std::vector<ad_scalar>& x,
                     std::vector<ad_scalar>& y) -> void {
    for (size_t j = 0; j < y.size(); ++j) y[j] = x[j];
  };
  t0 = clock::now();
  for (int i = 0; i < reps; ++i) {
    odelia::ode::vector_jacobian_product(tape, in, out_adjoint, nothing, in_adjoint);
  }
  const double t_tape =
    std::chrono::duration<double>(clock::now() - t0).count() / reps;

  return Rcpp::List::create(
    Rcpp::_["strategy_copy"] = t_strategy,
    Rcpp::_["environment_copy"] = t_environment,
    Rcpp::_["empty_tape_cycle"] = t_tape,
    Rcpp::_["sink"] = sink);
}

