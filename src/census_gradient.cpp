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
// accumulate in. Each name carries its species index, as "1.lma".
// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->r_patch().trait_adjoint_names();
}

// One row per census metric, one column per trait in census_trait_names_tf24()
// order. The active scalar lives inside this call and only doubles leave it.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_gradient_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_trait_gradient<plant::tf24_census>();
}

// The census's own reading of the traits at the state held. No sweep produces it,
// so it is the one route to the census a transpose check cannot touch, and it is
// easy to omit because it is a one-line calculation at the final state.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_direct_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_trait_direct<plant::tf24_census>();
}

// The same quantity differenced in plain double, which is what referees it.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_difference_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                             double rel) {
  return obj_->census_trait_difference<plant::tf24_census>(rel);
}

// The same gradient with the sweep stopped and resumed at each given recorded
// step, counted from one. Composition over steps is associative, so this must
// agree with census_trait_gradient_tf24 bit for bit.
// [[Rcpp::export]]
std::vector<std::vector<double>>
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
  return obj_->census_trait_gradient<plant::tf24_census>(at);
}

// How many backward ranges the last gradient swept. A requested split that fell
// on a segment boundary rather than inside one changes nothing, so this is what
// says a split actually cut.
// [[Rcpp::export]]
double census_adjoint_segments_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return static_cast<double>(obj_->adjoint_segments);
}

// [[Rcpp::export]]
std::vector<double> gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->gradient_control();
}
