#include <plant.h>

// The census exports stand_gradient() and stand_census() call.

// The census metrics and their trait gradient, for the one strategy that
// carries an active scalar. A metric is added to the strategy's own
// census_metrics(); nothing here names one, and the codomain is that list's
// length.

// [[Rcpp::export]]
std::vector<std::string> census_metric_names_tf24() {
  return plant::census_metric_names<plant::TF24_Strategy<double>>();
}

// [[Rcpp::export]]
std::vector<double> census_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census();
}

// [[Rcpp::export]]
std::vector<std::vector<double>>
census_state_adjoint_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_state_and_trait_rows().state.to_rows();
}

// The gradient's columns: each species' differentiable parameters in
// ad_parameters() order, species-major, which is the order the trait adjoints
// accumulate in. Each name carries its species index, as "1.lma".
// [[Rcpp::export]]
std::vector<std::string>
census_trait_names_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->r_patch().trait_adjoint_names();
}

// The gradient and its reading, as one R object. A caller cannot take the
// numbers without the statuses because they arrive in the same list: an exact
// zero and a refused row are both finite-looking, and only this says which.
Rcpp::List census_gradient_to_r(const plant::census_gradient& g) {
  const size_t n_metric = g.gradient.size();
  Rcpp::List gradient(n_metric), status(n_metric), refusal(n_metric);
  for (size_t m = 0; m < n_metric; ++m) {
    gradient[m] = Rcpp::wrap(g.gradient[m]);
    Rcpp::CharacterVector kinds(g.status[m].size());
    for (size_t p = 0; p < g.status[m].size(); ++p) {
      kinds[p] = plant::gradient_status::kind_name(g.status[m][p].kind);
    }
    status[m] = kinds;
    // Refusal is metric-level, so one description serves the row; a metric that
    // answered carries none rather than an empty one, so the two are not the
    // same object in R.
    const bool any_refused =
      !g.status[m].empty() &&
      g.status[m][0].kind == plant::gradient_status::Kind::refused;
    if (any_refused) {
      const plant::gradient_status& st = g.status[m][0];
      refusal[m] = Rcpp::List::create(
        Rcpp::_["reason"] = st.reason,
        Rcpp::_["species"] = st.species,
        Rcpp::_["node"] = st.node,
        Rcpp::_["step_first"] = st.step_first,
        Rcpp::_["step_last"] = st.step_last);
    } else {
      refusal[m] = R_NilValue;
    }
  }
  return Rcpp::List::create(Rcpp::_["gradient"] = gradient,
                            Rcpp::_["status"] = status,
                            Rcpp::_["refusal"] = refusal);
}

// One row per census metric, one column per trait in census_trait_names_tf24()
// order. The active scalar lives inside this call and only doubles leave it.
// [[Rcpp::export]]
Rcpp::List
census_trait_gradient_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_,
                           Rcpp::Nullable<Rcpp::CharacterVector> which_metrics = R_NilValue) {
  // Absent means every metric, which is what a caller that does not know the
  // census asks for. Naming a subset is what makes a single metric cost one, and
  // a name is what crosses: a position means a different metric as soon as the
  // model's list changes.
  std::vector<std::string> wanted;
  if (which_metrics.isNotNull()) {
    const Rcpp::CharacterVector v(which_metrics);
    for (int i = 0; i < v.size(); ++i) {
      wanted.push_back(Rcpp::as<std::string>(v[i]));
    }
  }
  return census_gradient_to_r(obj_->census_trait_gradient({}, wanted));
}

// [[Rcpp::export]]
std::vector<double> gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->gradient_control();
}
