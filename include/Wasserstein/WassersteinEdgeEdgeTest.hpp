// WassersteinEdgeEdgeTest.hpp

// include these from somewhere else

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <cmath>
#include <string>

using Real = float;

// exact predicates, but makes debugging and printing floating point in this routine hard
// using K = CGAL::Exact_predicates_exact_constructions_kernel;

// floating point predicates in the plane
using K = CGAL::Cartesian<Real>;

// real type underlying K
using K_real = K::RT;

using DT = CGAL::Delaunay_triangulation_2<K>;
using Point = DT::Point;
using Segment = CGAL::Segment_2<K>;

class WassersteinEdgeEdgeTest
{

public:

  // return the coordinates of the segment, as a nicely formatted string
  static
  std::string pretty_print( const Segment &e );

  // return the point i/(n-1) from e.source to e.target
  static
  Point frac(const Segment &e, int i, int n)
  {
    return e.source() + (i/ (double)(n-1)) * e.to_vector();
  }

  // compute W2 for two perpendicular segments, several different ways, to verify
  static
  double w2_perp(const Segment &e0, const Segment &e1);

  // compute W_p for two segments in general position, using brute force integration,
  // considering all permutations of the order of points
  // (more than 7 points might take a long time)
  // p is the order of W, typically p=1,2
  static
  double w_p(int p, const Segment &e0, const Segment &e1);

  // test the above for a bunch of different segments
  static
  void test();
};

// e.g.
/*
int main(int argc, char **argv) {
  
  // Test Wasserstein distance between two segments
  if (1)
  {
    WassersteinEdgeEdgeTest::test();
    return 0;
  }
}
*/
