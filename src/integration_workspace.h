// -*-c++-*-
#ifndef TREE_INTEGRATION_WORKSPACE_H_
#define TREE_INTEGRATION_WORKSPACE_H_

#include <list>
#include <vector>

namespace integration {

typedef std::vector< std::vector<double> > intervals_type;

namespace internal {

// These are used to simplify the book-keeping in the integrator, but
// could probably just be made part of the structure there?  However,
// keeping this out here means that the connections to the
// QUADPACK/GSL version is possibly a bit clearer.

class workspace {
public:
  class point {
  public:
    point(double a_, double b_, double area_, double error_)
      : a(a_), b(b_), area(area_), error(error_) {}
    double a;
    double b;
    double area;
    double error;
  };
  typedef std::list<point>::iterator points_iter;
  typedef std::list<point>::const_iterator points_const_iter;

  void clear();

  point worst_point() const;
  void update(point el1, point el2);

  double total_area() const;
  double total_error() const;

  void push_back(point el);
  intervals_type get_intervals() const;

private:
  void insert_forward(point el);
  void insert_backward(point el);
  void drop_worst();
  points_iter search_forward(double error);
  points_iter search_backward(double error);
  std::list<point> points;
};

intervals_type rescale_intervals(intervals_type x,
				 double min, double max);
}
}

#endif