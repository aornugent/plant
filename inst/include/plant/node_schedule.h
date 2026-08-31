// -*-c++-*-
#ifndef PLANT_PLANT_NODE_SCHEDULE_H_
#define PLANT_PLANT_NODE_SCHEDULE_H_

#include <RcppCommon.h> // SEXP
#include <plant/util.h>
#include <odelia/ode_interface.hpp>

// The "times" methods (set_times, times) refer to the *introduction*
// times.  As such, this really needs that at least one species has an
// introduction time of zero.  This is only enforced in that the model
// will just refuse to run if started with the incorrect time.

namespace plant {

// One introduction: the time, and the species that gain a node at it, ascending.
//
// A time naming several species is ONE introduction, because the run introduces
// them together and computes the environment once for the set. The type Patch
// takes is this species list, so the schedule hands over what introduce_nodes
// already wants.
struct introduction {
  double time;
  std::vector<size_t> species;
};

class NodeSchedule {
public:
  NodeSchedule(size_t n_species_);
  // Introductions, not species-times: a time naming three species is one.
  size_t size() const;
  size_t get_n_species() const;
  NodeSchedule expand(size_t n_extra, std::vector<double> times);
  void clear_times(size_t species_index);
  void set_times(const std::vector<double>& times_, size_t species_index);
  void set_times(const std::vector<std::vector<double> >& times);
  std::vector<double> times(size_t species_index) const;
  // Puts the walk back at the first introduction. The schedule itself is not
  // consumed by a run, so this restores a position and nothing else.
  void reset();

  // Where the walk is, where its interval ends, and moving past it. `next()`
  // stays valid across a pop, which only moves the position.
  const introduction& next() const;
  double time_end() const;
  void pop();
  size_t remaining() const;

  // The program a replay of the interval `(start, end)` takes: the state it
  // starts from, then the recorded steps strictly inside it. Empty where the
  // recording put no step in there, and empty on a schedule holding no recording.
  //
  // Boundaries are excluded because the run stops at them anyway and the
  // interval above starts from there -- so a step at one would be taken twice.
  std::vector<odelia::ode::instruction>
  program_within(double start, double end) const;

  double get_max_time() const;
  std::vector<std::vector<double> > get_times() const;
  bool using_ode_steps() const { return !ode_steps.empty(); }

  // * R interface:
  void r_clear_times(util::index species_index);
  std::vector<double> r_times(util::index species_index) const;
  void r_set_times(std::vector<double> times_, util::index species_index);
  void r_set_max_time(double x);
  std::vector<double> r_ode_times() const;
  std::vector<double> r_ode_step_sizes() const;
  // The two halves of one recording, installed together: apart, they can be
  // paired across different runs and nothing says so.
  void r_set_ode_steps(std::vector<double> times, std::vector<double> sizes);
  void r_clear_ode_steps();
  // Where the walk is, as a list: the time, the species from one, and the time
  // its interval ends. The species are one number each, counted the way R counts
  // them, where the class this replaces reported the same index twice under two
  // names.
  SEXP r_next_introduction() const;
  SEXP r_all_times() const;
  void r_set_all_times(SEXP x);
  NodeSchedule r_copy() const;

private:
  void insert(double time, size_t species_index);

  size_t n_species;
  // Ascending by time, one entry per distinct time, species ascending within an
  // entry -- which is the order Patch::introduced_at rebuilds off the same times,
  // so the two agree.
  std::vector<introduction> schedule;
  // How far the run has got. A position rather than a consumed copy of the
  // schedule, so what the run reads and what a caller set are one object.
  size_t at;
  double max_time;
  // The recorded ODE steps this schedule replays, when it holds any. Holding
  // them is what makes a schedule a replay, so no flag says so separately.
  std::vector<odelia::ode::instruction> ode_steps;
};

}

#endif
