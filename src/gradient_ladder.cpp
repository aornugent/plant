#include <plant.h>
#include <XAD/XAD.hpp>

// The reference the reverse sweep is checked against, and the two objects it is
// checked at.
//
// The reference is a forward-mode tangent of the same forward source. It is
// exact -- no step size, no truncation -- and it traverses the forward
// reductions while none of the transposes under test are on its path. One seed
// gives one exact Jacobian column, which is what makes forming a whole Jacobian
// affordable at one cohort.
//
// Two objects, because two different things are being refereed. The recorded
// cohort block relates one individual's states, the field's knot values and
// slopes, the soil potentials and the traits to its rates, its density rate and
// its per-layer draws; taking it forward and backward checks that the supplied
// leaf rows serve both modes. One right-hand-side evaluation relates the whole
// ODE state to its own rates through both reductions; taking it forward and
// backward checks the transposes, and there the reference and the object under
// test share no code.

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
  return Rcpp::List::create(Rcpp::_["state"] = lambda_y,
                            Rcpp::_["trait"] = patch.trait_adjoint,
                            Rcpp::_["block_recording_size"] =
                              static_cast<double>(patch.block_recording_size),
                            Rcpp::_["block_sweeps"] =
                              static_cast<double>(patch.block_sweeps));
}

// The trait columns' names, species-major, so a gradient's columns can be told
// apart when two species carry the same parameter.
// [[Rcpp::export]]
std::vector<std::string> ladder_trait_names_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  const patch_type& patch = *obj_;
  std::vector<std::string> out;
  for (size_t i = 0; i < patch.size(); ++i) {
    for (const std::string& n :
         patch.at_species(i).strategy_ptr()->ad_parameter_names()) {
      out.push_back(n);
    }
  }
  return out;
}

// How many nodes the block loop visits, which is what a node index runs over.
// [[Rcpp::export]]
int ladder_node_count_tf24(plant::RcppR6::RcppR6<plant::Patch<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return static_cast<int>(obj_->node_count());
}
