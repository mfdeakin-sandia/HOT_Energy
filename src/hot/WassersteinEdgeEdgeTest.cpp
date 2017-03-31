// WassersteinEdgeEdgeTest.cpp 

#include "WassersteinEdgeEdgeTest.hpp"

#include <algorithm>

#include <iostream>
#include <sstream>

std::string WassersteinEdgeEdgeTest::pretty_print( const Segment &e )
{
  std::stringstream ss;
  ss << "(" << e.source() << ")-(" << e.target() << ")";
  return ss.str();
}

Point WassersteinEdgeEdgeTest::frac(const Segment &e, int i, int n)
{
  return e.source() + (i/ (double)(n-1)) * e.to_vector();
}

double WassersteinEdgeEdgeTest::w2_perp(const Segment &e0, const Segment &e1)
{
  std::cout << "The Wasserstein2 distance between segements " << pretty_print(e0) << " and " << pretty_print(e1) << " is" << std::endl;

  const double c = sqrt(1/3.);

  // vectors between some pairs of the four points
  // ensure we match the closest ones
  auto ds = e0.source() - e1.source();
  auto dt = e0.target() - e1.target();
  
  // debug
  //  const double ds2 = CGAL::to_double( ds * ds );
  //  const double dt2 = CGAL::to_double( dt * dt );
  //  const double dst = CGAL::to_double( ds * dt );
  //  std::cout << "ds2:" << ds2 << " dt2:" << dt2 << " dst:" << dst << std::endl;
  
  const double W2_analytic = c * sqrt( CGAL::to_double( ds * ds + dt * dt + ds * dt ) );
  std::cout << "W2_analytic : " << W2_analytic << std::endl;
  
  // compare vs. matching the endpoints the other way
  // this only makes a difference for oblique (non-perpendicular, non-parallel) edges.
  if (1)
  {
    auto dsR = e0.source() - e1.target();
    auto dtR = e0.target() - e1.source();
    const double W2_analyticR = c * sqrt( CGAL::to_double( dsR * dsR + dtR * dtR + dsR * dtR ) );
    std::cout << "W2_anal_Rev : " << W2_analyticR << std::endl;
  }
  
  
  // verify by doing three dot products manually
  if (1)
  {
    K_real ddss(0.), ddtt(0), ddst(0);
    for (int i = 0; i < ds.dimension(); ++i)
    {
      ddss  += ds.cartesian(i) * ds.cartesian(i);
      ddtt  += dt.cartesian(i) * dt.cartesian(i);
      ddst  += ds.cartesian(i) * dt.cartesian(i);
    }
    const double W2_analytic_manual = c * sqrt( CGAL::to_double( ddss + ddtt + ddst ) );
    std::cout << "W2_manual   : " << W2_analytic_manual << std::endl;
  }

  
  // verification against numeric integration
  if (1)
  {
    // numeric
    //    for (int n = 10; n < 123456; n *= 10)
    for (int n = 10000; n < 12345; n *= 10)
    {
      K_real d_2(0.);
      for (int i = 0; i < n; ++i )
      {
        Point p0 = frac( e0, i, n );
        Point p1 = frac( e1, i, n );
        d_2 += CGAL::squared_distance(p0,p1);
      }
      const double W2_numeric = sqrt( CGAL::to_double(d_2) / (double) n );
      
      std::cout << "W2_numeric  : " << W2_numeric << std::endl;
      std::cout << "W2_numeric^2 / W2_analytic^2 = " << W2_numeric * W2_numeric / (W2_analytic * W2_analytic) << std::endl;
      std::cout << "W2_analytic^2 / W2_numeric^2 = " << (W2_analytic * W2_analytic) / (W2_numeric * W2_numeric) << std::endl;
    }
  }
  
  return W2_analytic;
}

double WassersteinEdgeEdgeTest::w_p(int p, const Segment &e0, const Segment &e1)
{
  // std::cout << "Computing the Wasserstein distance between segements " << pretty_print(e0) << " and " << pretty_print(e1) << std::endl;
  
  // dumb brute force way
  // each integration point of one segment is mapped to an integration point in another

  // n = number of integration points,
  // if using exact,   n <= 7 for fast, n = 8 involves 15 second wait
  // if using Cartesian, n=8 is instant, n = 10 takes 30 seconds,
  const int n=6;
  assert( n > 1 );
  assert( p >= 0 );

  // data for permutations
  // numbers from 0 to n-1
  std::vector<int> ind(n);
  std::iota(ind.begin(), ind.end(), 0);
  std::vector<int> best_ind(ind);
  std::vector<int> worst_ind(ind);

  // consider all possible pairings of integration points
  // go through all the points of e0 in order
  // go through all the points of e1 in all possible permutations

  double min_dp_double = std::numeric_limits<double>::max();
  double worst_dp_double = 0.;
  K_real min_d2( min_dp_double );
  K_real worst_d2( worst_dp_double );

  // for all permutations
  do
  {
    
    K_real d2(0.);
    double dp_double(0.);
    
    // all points of e0
    for (auto i = 0; i < n; ++i)
    {
      Point p0 = frac( e0, i, n );
      Point p1 = frac( e1, ind[i], n );
    
      auto squared_d = CGAL::squared_distance( p0, p1 );  // W2
      if (p==2)
        d2 += squared_d;
      else
      {
        dp_double += pow( CGAL::to_double( squared_d ), (int) (p / 2) );  // W_p
      }

//      std::cout << "i: " << i << " p:" << p0 << " to j:" << ind[i] << " q:" << p1 << " d^2:" << d2 << " == " << d2_float << std::endl;
    }
//    std::cout << "d2:" << d2 << std::endl;
//     std::cout << "  " << d2;
    
    // d2 = sum of distances^2 for this permutation
    // best so far?
    zzyk
    if ( (p==2 ? d2 < min_d2 : dp_double < min_dp_double ) )
    {
      best_ind = ind;
      min_d2 = d2;
//      std::cout << "new min :" << min_d2 << std::endl;
    }
    if ( (p==2 ? d2 > worst_d2 : dp_double > worst_dp_double ) )
    {
      worst_ind = ind;
      worst_d2 = d2;
//      std::cout << "new worst :" << worst_d2 << std::endl;
    }

  } while ( std::next_permutation(ind.begin(), ind.end()) );
  // if permutation speed becomes important, try http://howardhinnant.github.io/combinations.html
  
  auto min_d = sqrt( CGAL::to_double(min_d2 / n) ); // W2
  auto worst_d = sqrt( CGAL::to_double(worst_d2 / n) ); // W2
  // auto min_d = CGAL::to_double(min_d2) / (double) n; // W1
  // auto worst_d = CGAL::to_double(worst_d2) / (double) n; // W1
  
  std::cout << "\nThe Wasserstein2 distance between segements " << pretty_print(e0) << " and " << pretty_print(e1) << " was " << min_d << std::endl;
  std::cout << "From permutation  :";
  for ( auto p : best_ind )
    std::cout << " " << p;
  std::cout << "\nWorst permutation :";
  for ( auto p : worst_ind )
    std::cout << " " << p;
  std::cout << " gave distance " << worst_d << ". ";
  
  std::cout << std::endl << std::endl;

  return min_d;
}

void wasserstein_edge_edge_test()
{
  // create two test edges
  Point p0(0,0);
  Point p1(2,0);
  Point p2(1,0);
  Point p3(1,2);
  
  // use P10--P11 for both, so distance should be zero
  Segment e0(p0,p1);
  Segment e1(p0,p1);

  // compute distance
//  wasserstein2_edge_edge( Segment(p0,p1), Segment(p0,p1) );
//  wasserstein2_edge_edge( Segment(p0,p1), Segment(p1,p0) );
  
  // perp edges at midpoint
  // of uniform length
  // appears that W2 is invariant to permutation of matching up integration points!
  // W1 order matters, best is a large spread, match close points together and match far points together
  // W1 best is interlaced, alternating points from each side, *not* 0-5 matched with 0-5
  if (1)
  {
    Segment s01( Point(0,0), Point(1,0) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5, 0.0), Point(0.5,0.0)) );
    
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5, 1.0), Point(0.5,2.0)) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5, 0.5), Point(0.5,1.5)) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5, 0.0), Point(0.5,1.0)) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5,-0.2), Point(0.5,0.8)) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5,-0.4), Point(0.5,0.6)) );
    wasserstein2_edge_edge_perp( s01, Segment( Point(0.5,-0.5), Point(0.5,0.5)) );
  }
  if (1)
  {
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1,1), Point(1,3)) );      // 2.22  W2
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1,0), Point(1,2)) );      // 1.39
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1,-0.3), Point(1,1.7)) ); // 1.19
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1,-0.5), Point(1,1.5)) ); // 1.09
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1,-1), Point(1,1)) );     // 0.966
  }
  
  // perp edges, not at midpoint
  // appears that W2 is invariant to permutation of matching up integration points!
  // W1 is as in the above case
  if (1)
  {
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,1), Point(0,3)) );      // 2.43 W2
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,0), Point(0,2)) );      // 1.71
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-0.3), Point(0,1.7)) ); // 1.55
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-0.5), Point(0,1.5)) ); // 1.47
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-1), Point(0,1)) );     // 1.39
  }
  
  // oblique
  if (1)
  {
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(-1,1), Point(1,1)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(-1,1), Point(1,2)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(-1,1), Point(1,3)) );

    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,1)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,2)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,3)) );

    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.8,0.3)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.8,0.9)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.5,1.8)) );
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(1.5,-0.3), Point(0.2,1.8)) );
  }
}

// e.g.
/*
int main(int argc, char **argv) {
  
  // Test Wasserstein distance between two segments
  if (1)
  {
    wasserstein_edge_edge_test();
    return 0;
  }
}
*/
