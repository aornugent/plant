#include <plant.h>

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
  return obj_->census_state_adjoint().to_rows();
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
namespace {
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
}  // namespace

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

// The census's own reading of the traits at the state held. No sweep produces it,
// so it is the one route to the census a transpose check cannot touch, and it is
// easy to omit because it is a one-line calculation at the final state.
// [[Rcpp::export]]
std::vector<std::vector<double>>
census_trait_direct_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->census_trait_direct().to_rows();
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

// [[Rcpp::export]]
std::vector<double> gradient_control_tf24(plant::RcppR6::RcppR6<plant::SCM<plant::TF24_Strategy<double>, plant::TF24_Environment<double> > > obj_) {
  return obj_->gradient_control();
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
