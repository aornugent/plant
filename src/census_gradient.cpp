#include <plant.h>

// The census metrics and their trait gradient, for the one strategy that
// carries an active scalar. A metric is added in species.h; nothing here names
// one, so the codomain follows the tuple.

// [[Rcpp::export]]
std::vector<std::string> census_metric_names_tf24() {
  std::vector<std::string> ret;
  std::apply([&](auto... psi) -> void { (ret.push_back(psi.name()), ...); },
             plant::tf24_census{});
  return ret;
}

// [[Rcpp::export]]
std::vector<double> census_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census<plant::tf24_census>();
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
census_state_adjoint_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_state_adjoint<plant::tf24_census>();
}

// The gradient's columns: each species' differentiable parameters in
// ad_parameters() order, species-major, which is the order the trait adjoints
// accumulate in.
// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  std::vector<std::string> ret;
  const plant::Patch<plant::TF24_Strategy<double>,
                     plant::TF24_Environment<double>>& patch = obj_->r_patch();
  for (size_t i = 0; i < patch.size(); ++i) {
    for (const std::string& n :
         patch.at_species(i).strategy_ptr()->ad_parameter_names()) {
      ret.push_back(n);
    }
  }
  return ret;
}

// [[Rcpp::export]]
std::vector<double> gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->gradient_control();
}
