// WassersteinEdgeEdge.hpp

// include these from somewhere else

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <cmath>

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

// This works for a 2d mesh embedded in the plane, or a 3d mesh in 3d space, etc.

// For 2d mesh embedded in 3d, we will have two dual half-edges, 
// one in the plane of each triangle of the edge.
// Compute the distance (which could be negative) for each half-edge separately and add, 
// as in the paper
// Technical note
// Delaunay Hodge star
// Anil N. Hirani a,∗, Kaushik Kalyanaraman a, Evan B. VanderZee b

class WassersteinEdgeEdge
{
public:
  // W2 (Wasserstein p=2) distance between two segments
  // edge e0 = (e0s, e0t) // two points of edge
  // edge e1 = (e1s, e1t)
  // vector ds = e0s - e1s
  // vector dt = e0t - e1t
  // distance^2 = (ds * ds + dt * dt + ds * dt) / 3 // * denotes dot product
  // if e0 is oblique (not-perpendicular and not-parallel) to e1, then 
  //   recompute with reversing the direction of e1 and e1 
  //   return the smaller of the two distances

  static double w2_perp_squared(const Segment &e0, const Segment &e1)
  {
    // vectors between some pairs of the four points
    // for perpendicular edges, including perpendicular bisectors, it doesn't matter which pair we match
    const auto ds = e0.source() - e1.source();
    const auto dt = e0.target() - e1.target();
    // constant of integration
    const double c = 1./3.;
    const double W2_analytic_squared = c * CGAL::to_double( ds * ds + dt * dt + ds * dt );
    return W2_analytic_squared;
  }  
  static double w2_perp(const Segment &e0, const Segment &e1)
  {
    return sqrt( w2_perp_squared );
  }


  static double w2_oblique_squared(const Segment &e0, const Segment &e1)
  {
    // try both directions for e1 ensures we match the orientation giving the shortest distance
    // since lines are straight, for W2 there is no need to try anything else
    //                           for W1 some mixed order may be smaller
    const double W2_forward_2 = w2_perp_squared( e0, e1 );
    const Segment e1_reverse( e1.target(), e1.source() );
    const double W2_reverse_2 = W2_perp_squared( e0, e1_reverse );
    return std::min( W2_forward, W2_reverse );
  }
  static double w2_oblique(const Segment &e0, const Segment &e1)
  {
    return sqrt( w2_oblique_squared );    
  }

  // W1_perp (bisector)
  // W1_oblique <- skip. This is complicated, because of the matching order.

};
