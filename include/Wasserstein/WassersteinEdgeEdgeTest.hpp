// WassersteinEdgeEdgeTest.hpp

#ifndef WASSERSTEINEDGEEDGETEST_HPP
#define WASSERSTEINEDGEEDGETEST_HPP

#include "WassersteinKernel.hpp"

class WassersteinEdgeEdgeTest
{

public:

  // return the coordinates of the segment, as a nicely formatted string
  static
  std::string pretty_print( const Segment &e );

  static
  Point frac(const Segment &e, double f)
  {
    return e.source() + f * e.to_vector();
  }
  // return the point i/(n-1) from e.source to e.target
  static
  Point frac(const Segment &e, int i, int n)
  {
    return frac( e, (i/ (double)(n-1)) );
  }

  // compute W1 "efficiently"
  // for an edge and its dual, 
  // i.e. for two perpendicular segments, where the affine hull of one passes through the midpoint of the other,
  // Compute several different ways for verification.
  // The numerical verification uses the permutation, so isn't very accurate and
  // is limited to only a few integration points, but seems to converge to what I'm calculating.
  // The numerical integration always gives a larger answer than the fast computation,
  // which is consistent with the matching being poor for a coarse set of points.
  static 
  double w1_edge_dual(const Segment &e0, Segment e1);

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

#endif // WASSERSTEINEDGEEDGETEST_HPP
