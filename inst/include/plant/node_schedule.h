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

// This could be done via a list object, but I think this is OK for
// now.  The main reason for keeping this as a separate class is it
// only makes sense to have a nontrivial constructor, and that's not
// yet supported for RcppR6 lists.
class NodeScheduleEvent {
public:
  NodeScheduleEvent(double introduction, size_t species_index_)
    : species_index(species_index_) {
    times.push_back(introduction);
  }
  size_t species_index_raw() const {
    return species_index;
  }
  double time_introduction() const {
    return times.front();
  }
  double time_end() const {
    return times.back();
  }

  size_t species_index;
  std::vector<double> times;
  // The ODE steps a recorded run took inside this interval, each carrying the
  // time it reached and the size that reached it. Empty unless the schedule
  // holds a recording.
  std::vector<odelia::ode::recorded_step> steps;
};

class NodeSchedule {
public:
  typedef NodeScheduleEvent Event;
  NodeSchedule(size_t n_species_);
  size_t size() const;
  size_t get_n_species() const;
  NodeSchedule expand(size_t n_extra, std::vector<double> times);
  void clear_times(size_t species_index);
  void set_times(const std::vector<double>& times_, size_t species_index);
  void set_times(const std::vector<std::vector<double> >& times);
  std::vector<double> times(size_t species_index) const;
  void reset();
  void pop();
  Event next_event() const;
  size_t remaining() const;

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
  SEXP r_all_times() const;
  void r_set_all_times(SEXP x);
  NodeSchedule r_copy() const;

private:
  typedef std::list<Event>::iterator events_iterator;
  typedef std::list<Event>::const_iterator events_const_iterator;

  events_iterator add_time(double times, size_t species_index,
                           events_iterator it);
  void distribute_ode_steps();

  size_t n_species;
  std::list<Event> events;
  std::list<Event> queue;
  double max_time;
  // The run this schedule replays, if it holds one. A schedule with a recording
  // replays it; there is no third state and so no flag. What used to be two
  // vectors and a bool was one recording split three ways, held together by a
  // rule about which setter cleared which.
  std::vector<odelia::ode::recorded_step> ode_steps;
};

}

#endif
