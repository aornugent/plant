#include <plant/node_schedule.h>
#include <plant/parameters.h>
#include <plant/util.h>
#include <Rcpp.h>
#include <cmath> // log2, exp2
#include <limits> // std::numeric_limits

namespace plant {

NodeSchedule::NodeSchedule(size_t n_species_)
  : n_species(n_species_),
    max_time(std::numeric_limits<double>::infinity()) {
  reset();
}

size_t NodeSchedule::size() const {
  return events.size();
}

size_t NodeSchedule::get_n_species() const {
  return n_species;
}

NodeSchedule NodeSchedule::expand(size_t n_extra,
                                  std::vector<double> times) {
  NodeSchedule ret = *this;
  ret.n_species += n_extra;
  for (size_t i = n_species; i < ret.n_species; ++i) {
    ret.set_times(times, i);
  }
  return ret;
}

void NodeSchedule::clear_times(size_t species_index) {
  events_iterator e = events.begin();
  while (e != events.end()) {
    if (e->species_index == species_index) {
      e = events.erase(e);
    } else {
      ++e;
    }
  }
  reset();
}

double NodeSchedule::get_max_time() const {
  return max_time;
}

std::vector<std::vector<double> > NodeSchedule::get_times() const {
  std::vector<std::vector<double> > ret;
  for (size_t i = 0; i < n_species; ++i) {
    ret.push_back(times(i));
  }
  return ret;
}

void NodeSchedule::set_times(const std::vector<double>& times_,
                               size_t species_index) {
  clear_times(species_index);
  events_iterator e = events.begin();
  for (std::vector<double>::const_iterator t = times_.begin();
       t != times_.end(); ++t) {
    e = add_time(*t, species_index, e);
  }
  reset();
}

void NodeSchedule::set_times(const std::vector<std::vector<double> >& times_) {
  util::check_length(times_.size(), n_species);
  for (size_t i = 0; i < n_species; ++i) {
    set_times(times_[i], i);
  }
}

std::vector<double> NodeSchedule::times(size_t species_index) const {
  std::vector<double> ret;
  for (events_const_iterator e = events.begin(); e != events.end(); ++e) {
    if (e->species_index == species_index) {
      ret.push_back(e->time_introduction());
    }
  }
  return ret;
}

void NodeSchedule::reset() {
  queue = events;

  // NOTE: this is updating elements in *queue*, not in events; the
  // events only ever have a single element in times.
  events_iterator e = queue.begin();
  while (e != queue.end()) {
    events_iterator e_next = e;
    ++e_next;
    e->times.push_back(e_next == queue.end() ?
                       max_time : e_next->time_introduction());
    e = e_next;
  }

  if (using_ode_steps()) {
    distribute_ode_steps();
  }
}

// NOTE: Using vector here is inefficient, but we need a vector in the
// end.  This would not matter if we added to event and didn't do this
// each reset, but I don't imagine that this is a big cost in
// practical cases.
//
// NOTE: It should be the case that there are exactly two times in
// e->times on entry, which means that we could always just construct
// a new vector with the start time, the times in 'extra' and the end
// time, purely by pushing on the end.  Ignoring this, it does mean
// that the inefficiency of using vector (over list) is confined to a
// copy of a single element, which won't be as bad as using an
// underlying list and converting to vector when passing back to the
// ode solver.
void NodeSchedule::distribute_ode_steps() {
  // Each recorded step carries the time it reached, so an interval's steps need
  // no second vector to stay aligned with anything.
  events_iterator e = queue.begin();
  size_t i = 0;
  while (e != queue.end()) {
    std::vector<odelia::ode::recorded_step> inside;
    while (i < ode_steps.size() && ode_steps[i].time < e->time_end()) {
      // Excludes times that exactly match an interval boundary; the run stops
      // there anyway, and the interval above it starts from there.
      if (!util::identical(ode_steps[i].time, e->time_introduction()) &&
          !util::identical(ode_steps[i].time, e->time_end())) {
        inside.push_back(ode_steps[i]);
      }
      ++i;
    }
    if (!inside.empty()) {
      e->steps.assign(1, {e->time_introduction(),
                          std::numeric_limits<double>::quiet_NaN()});
      e->steps.insert(e->steps.end(), inside.begin(), inside.end());
      std::vector<double> extra;
      extra.reserve(inside.size());
      for (const odelia::ode::recorded_step& r : inside) {
        extra.push_back(r.time);
      }
      std::vector<double>::iterator at = ++e->times.begin();
      e->times.insert(at, extra.begin(), extra.end());
    }
    ++e;
  }
}

void NodeSchedule::pop() {
  if (queue.empty()) {
    Rcpp::stop("Attempt to pop empty queue");
  }
  queue.pop_front();
}

NodeScheduleEvent NodeSchedule::next_event() const {
  if (queue.empty()) {
    Rcpp::stop("All events completed");
  }
  return queue.front();
}

size_t NodeSchedule::remaining() const {
  return queue.size();
}

// * R interface
void NodeSchedule::r_clear_times(util::index species_index) {
  clear_times(species_index.check_bounds(n_species));
}

void NodeSchedule::r_set_times(std::vector<double> times_,
                                 util::index species_index) {
  if (!util::is_sorted(times_.begin(), times_.end())) {
    Rcpp::stop("Times must be sorted (increasing)");
  }
  if (times_.size() == 0) {
    Rcpp::stop("Need at least one time");
  }
  if (times_.front() < 0) {
    Rcpp::stop("First time must nonnegative");
  }
  if (times_.back() > max_time) {
    Rcpp::stop("Times cannot be greater than max_time");
  }
  set_times(times_, species_index.check_bounds(n_species));
}

std::vector<double> NodeSchedule::r_times(util::index species_index) const {
  return times(species_index.check_bounds(n_species));
}

void NodeSchedule::r_set_max_time(double x) {
  if (x < 0) {
    Rcpp::stop("max_time must be nonnegative");
  }
  if (x < events.back().time_introduction()) {
    Rcpp::stop("max_time must be at least the final scheduled time");
  }
  max_time = x;
  reset();
}

std::vector<double> NodeSchedule::r_ode_times() const {
  std::vector<double> ret;
  ret.reserve(ode_steps.size());
  for (const odelia::ode::recorded_step& r : ode_steps) {
    ret.push_back(r.time);
  }
  return ret;
}

std::vector<double> NodeSchedule::r_ode_step_sizes() const {
  std::vector<double> ret;
  ret.reserve(ode_steps.size());
  for (const odelia::ode::recorded_step& r : ode_steps) {
    ret.push_back(r.step_size);
  }
  return ret;
}

// The ODE schedule this run takes: times to stop at, and the sizes that reached
// them where a recorded run is being replayed. Both at once, because a size
// belongs to the time it was recorded with and there is no version of this that
// takes one and keeps the other.
//
// Sizes may be omitted, and that is a grid a caller chose rather than a run a
// caller recorded: the solver steps TO each time instead of by a recorded size.
// A schedule holding either replays it -- there is no flag to turn that on, so to
// integrate freely instead, hand over no schedule.
void NodeSchedule::r_set_ode_steps(std::vector<double> times,
                                   std::vector<double> sizes) {
  if (times.empty()) {
    r_clear_ode_steps();
    return;
  }
  if (times.size() < 2) {
    Rcpp::stop("Need at least two times");
  }
  if (sizes.empty()) {
    // A grid: every time is stepped to, so no size is known anywhere.
    sizes.assign(times.size(), std::numeric_limits<double>::quiet_NaN());
  }
  if (sizes.size() != times.size()) {
    Rcpp::stop("ode_step_sizes must be the same length as ode_times");
  }
  if (!util::identical(times.front(), 0.0)) {
    Rcpp::stop("First time must be exactly zero");
  }
  if (util::is_finite(max_time) && !util::identical(times.back(), max_time)) {
    Rcpp::stop("Last time must be exactly max_time");
  }
  if (!util::is_sorted(times.begin(), times.end())) {
    Rcpp::stop("ode_times must be sorted");
  }
  if (!std::isnan(sizes.front())) {
    Rcpp::stop("First step size must be NaN, the recorded start that no step "
               "reached");
  }
  ode_steps.clear();
  ode_steps.reserve(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    ode_steps.push_back({times[i], sizes[i]});
  }
  if (!util::is_finite(max_time)) {
    max_time = ode_steps.back().time;
  }
  reset();
}

void NodeSchedule::r_clear_ode_steps() {
  ode_steps.clear();
  // This will pull the recorded steps out of the events if they were set.
  reset();
}


SEXP NodeSchedule::r_all_times() const {
  return Rcpp::wrap(get_times());
}

void NodeSchedule::r_set_all_times(SEXP rx) {
  Rcpp::List x(Rcpp::as<Rcpp::List>(rx));
  // Ensure that we can get all the times out:
  std::vector< std::vector<double> > new_times;
  for (Rcpp::List::iterator el = x.begin(); el != x.end(); ++el) {
    new_times.push_back(Rcpp::as<std::vector<double> >(*el));
  }
  // set_all_times(new_times);
  util::check_length(new_times.size(), n_species);
  for (size_t i = 0; i < n_species; ++i) {
    set_times(new_times[i], i);
  }
}

NodeSchedule NodeSchedule::r_copy() const {
  return *this;
}

// * Private methods
NodeSchedule::events_iterator
NodeSchedule::add_time(double time, size_t species_index,
                         events_iterator it) {
  Event e(time, species_index);
  it = events.begin();
  while (it != events.end() && time > it->time_introduction()) {
    ++it;
  }
  it = events.insert(it, e);
  return it;
}

}
