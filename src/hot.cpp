// hot.cpp

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <iostream>

#include <polynomial.hpp>

using RNG = std::mt19937_64;

using Real = float;

constexpr const int dims = 2;

// exact predicates, but makes debugging and printing floating point in this routine hard
// using K = CGAL::Exact_predicates_exact_constructions_kernel;

// floating point predicates in the plane
using K = CGAL::Cartesian<Real>;

// real type underlying K
using K_real = K::RT;

using DT = CGAL::Delaunay_triangulation_2<K>;
using Face = DT::Face;
using Point = DT::Point;

using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

constexpr const int tri_verts = 3;
constexpr const int tri_edges = 3;

constexpr const float min_pos = 10.0;
constexpr const float max_pos = 20.0;

void generate_rand_dt(int num_points, DT &dt) {
  std::random_device rd;
  RNG engine(rd());
  std::uniform_real_distribution<Real> genPos(min_pos,
                                              max_pos);
  for(int i = 0; i < num_points; i++) {
    Point pt(genPos(engine), genPos(engine));
    dt.insert(pt);
  }
}

K::RT triangle_area(const Face &face) {
  return Triangle(face.vertex(0)->point(),
                  face.vertex(1)->point(),
                  face.vertex(2)->point())
      .area();
}

/* Computes the centroid of the triangle */
Point triangle_centroid(const Face &face) {
  return CGAL::centroid(face.vertex(0)->point(),
                        face.vertex(1)->point(),
                        face.vertex(2)->point());
}

/* Computes 7 triangles which can be used for the piecewise
 * integral of the Wasserstein distance with odd k */
std::array<Triangle, 7> odd_integral_bounds(
    const Face &face, const Point &centroid) {
  const Line vertical(centroid, centroid + Vector(0, 1));
  const Line horizontal(centroid, centroid + Vector(1, 0));
  Segment edges[tri_edges];
  for(int i = 0; i < tri_edges; i++) {
    edges[i] = Segment(face.vertex(face.ccw(i))->point(),
                       face.vertex(face.cw(i))->point());
  }
  for(int i = 0; i < tri_edges; i++) {
    auto int_vert = CGAL::intersection(vertical, edges[i]);
    // boost::optional - We need to verify that there is an
    // intersection
    if(int_vert.is_initialized()) {
      // Any non-degenerate triangle will cause this
      // intersection to be a point.
      Point p = boost::get<Point>(int_vert.get());
    }

    auto int_horiz =
        CGAL::intersection(horizontal, edges[i]);
    // boost::optional - We need to verify that there is an
    // intersection
    if(int_horiz.is_initialized()) {
      // Any non-degenerate triangle will cause this
      // intersection to be a point.
      Point p = boost::get<Point>(int_vert.get());
    }
  }
}

void order_points(std::array<Point, tri_verts> &verts) {
  if(verts[0][0] > verts[1][0]) {
    std::swap(verts[0], verts[1]);
  }
  if(verts[1][0] > verts[2][0]) {
    std::swap(verts[2], verts[1]);
    if(verts[0][0] > verts[1][0]) {
      std::swap(verts[0], verts[1]);
    }
  }
}

/* Computes 1 or 2 triangles which can be used for the
 * piecewise integral of the Wasserstein distance with even
 * k */
boost::variant<std::array<Triangle, 2>,
               std::array<Triangle, 1> >
even_integral_bounds(const Face &face) {
  std::array<Point, tri_verts> verts = {
      face.vertex(0)->point(), face.vertex(1)->point(),
      face.vertex(2)->point()};
  order_points(verts);
  if(verts[0][0] == verts[1][0] ||
     verts[1][0] == verts[2][0]) {
    // Return a single triangle in this case
    std::array<Triangle, 1> bounds;
    bounds[0] = Triangle(verts[0], verts[1], verts[2]);
    return boost::variant<std::array<Triangle, 2>,
                          std::array<Triangle, 1> >(bounds);
  } else {
    const Line vertical(verts[1], verts[1] + Vector(0, 1));
    const Line base(verts[0], verts[2]);
    // This intersection will exist for all non-degenerate
    // triangles
    auto int_vert = CGAL::intersection(vertical, base);
    assert(int_vert.is_initialized());

    std::array<Triangle, 2> bounds;
    bounds[0] = Triangle(verts[0], verts[1],
                         boost::get<Point>(int_vert.get()));
    bounds[1] = Triangle(verts[2], verts[1],
                         boost::get<Point>(int_vert.get()));
    return boost::variant<std::array<Triangle, 2>,
                          std::array<Triangle, 1> >(bounds);
  }
}

template <typename T, int k>
class triangle_w_helper;

template <typename T>
class triangle_w_helper<T, 2> {
 public:
  K::RT operator()(T &bounds) {
    K::RT integral = 0;
    for(auto area : bounds) {
      /* Given a signed width of w, an initial x of x_0, two
       * bounding segments with y = si x + offi, and a
       * centroid located at (x_c, y_c), we need to compute:
       * w**4*(-s1**3 - 3*s1 + s2**3 + 3*s2)/12 +
       * w**3*(-off1*s1**2 - off1 + off2*s2**2 + off2 +
       * s1**2*y_c + 2*s1*x_c - s2**2*y_c - 2*s2*x_c)/3 +
       * w**2*(-off1**2*s1 + 2*off1*s1*y_c + 2*off1*x_c +
       * off2**2*s2 - 2*off2*s2*y_c - 2*off2*x_c - s1*x_c**2
       * - s1*y_c**2 + s2*x_c**2 + s2*y_c**2)/2 - w*(off1**3
       * - 3*off1**2*y_c + 3*off1*x_c**2 + 3*off1*y_c**2 -
       * off2**3 + 3*off2**2*y_c - 3*off2*x_c**2 -
       * 3*off2*y_c**2)/3 + x_0**4*(s1**3 + 3*s1 - s2**3 -
       * 3*s2)/12 + x_0**3*(off1*s1**2 + off1 - off2*s2**2 -
       * off2 - s1**2*y_c - 2*s1*x_c + s2**2*y_c +
       * 2*s2*x_c)/3 + x_0**2*(off1**2*s1 - 2*off1*s1*y_c -
       * 2*off1*x_c - off2**2*s2 + 2*off2*s2*y_c +
       * 2*off2*x_c + s1*x_c**2 + s1*y_c**2 - s2*x_c**2 -
       * s2*y_c**2)/2 + x_0*(off1**3 - 3*off1**2*y_c +
       * 3*off1*x_c**2 + 3*off1*y_c**2 - off2**3 +
       * 3*off2**2*y_c - 3*off2*x_c**2 - 3*off2*y_c**2)/3
       *
       * Start by computing the lines used for boundaries
       * and the width.
       */
      std::array<Point, tri_verts> verts = {
          area.vertex(0), area.vertex(1), area.vertex(2)};
      order_points(verts);

      /* Find the index of the non-vertical point
                         * It's guaranteed to be either the
       * leftmost or
                         * rightmost point
                         */
      const int w_idx =
          (verts[0][0] != verts[1][0]) ? 0 : 2;
      // Compute the signed width
      K::RT width = verts[w_idx][0] -
                    verts[(w_idx + 1) % tri_verts][0];
      /* Compute the bounding lines
       * x = ax t + bx -> t = (x - bx) /
       * ax
       * y = ay t + by = (ay / ax) x +
       * (-ay bx / ax + by)
       *
       * bx = verts[w_idx][0], by = verts[w_idx][1]
       */
      Vector dir_1 =
          Line(verts[w_idx], verts[(w_idx + 1) % tri_verts])
              .to_vector();
      Vector dir_2 =
          Line(verts[w_idx], verts[(w_idx + 2) % tri_verts])
              .to_vector();
      K::RT slope_1 = dir_1[1] / dir_1[0];
      K::RT slope_2 = dir_2[1] / dir_2[0];
      if((w_idx == 0 && slope_1 < slope_2) ||
         (w_idx == 2 && slope_1 > slope_2)) {
        std::swap(slope_1, slope_2);
      }
      Numerical::Polynomial<K::RT, 1, 1> bound_1;
      bound_1.coeff(1) = slope_1;
      bound_1.coeff(0) =
          -verts[w_idx][0] * slope_1 + verts[w_idx][1];
      Numerical::Polynomial<K::RT, 1, 1> bound_2;
      bound_2.coeff(1) = slope_2;
      bound_2.coeff(0) =
          -verts[w_idx][0] * slope_2 + verts[w_idx][1];

      Numerical::Polynomial<K::RT, 2, 2> initial(
          (Tags::Zero_Tag()));
      initial.coeff(2, 0) = 1;
      initial.coeff(0, 2) = 1;
      auto y_int = initial.integrate(1);
      Numerical::Polynomial<K::RT, 3, 1> upper =
          y_int.var_sub(1, 0, bound_1.coeff(1)) +
          y_int.slice(1, bound_1.coeff(0));
      Numerical::Polynomial<K::RT, 3, 1> lower =
          y_int.var_sub(1, 0, bound_2.coeff(1)) +
          y_int.slice(1, bound_2.coeff(0));
      auto y_bounded = upper + -lower;
      auto x_int = y_bounded.integrate(0);

      auto left = x_int.slice(0, verts[w_idx][0]).coeff(0);
      auto right =
          x_int.slice(0, verts[(w_idx + 1) % tri_verts][0])
              .coeff(0);
      integral += left + -right;
    }
    return integral;
  }
};

/* Computes the k Wasserstein distance of the triangular
 * face to it's centroid */
template <int k>
K::RT triangle_w(const Face &face);

template <>
K::RT triangle_w<2>(const Face &face) {
  Point centroid = triangle_centroid(face);
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = even_integral_bounds(face);
  // Set to NaN until implemented
  Real distance = 0.0 / 0.0;
  if(bounds.which() == 0) {
    // std::array<Triangle, 2>
    auto area =
        boost::get<std::array<Triangle, 2> >(bounds);
    return triangle_w_helper<std::array<Triangle, 2>, 2>()(
        area);
  } else {
    // std::array<Triangle, 1>
    auto area =
        boost::get<std::array<Triangle, 1> >(bounds);
    return triangle_w_helper<std::array<Triangle, 1>, 2>()(
        area);
  }
}

/* See "HOT: Hodge-Optimized Triangulations" for details on
 * the energy functional */
K::RT hot_energy(const DT &dt) {
  K::RT energy = 0;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end(); face_itr++) {
    K::RT area = triangle_area(*face_itr);
    K::RT wasserstein = triangle_w<2>(*face_itr);
    energy += wasserstein * area;
  }
  return energy;
}

std::string pretty_print( const Segment &e )
{
  std::stringstream ss;
  ss << "(" << e.source() << ")-(" << e.target() << ")";
  return ss.str();
}

Point frac(const Segment &e, int i, int n)
{
  return e.source() + (i/ (double)(n-1)) * e.to_vector();
}

double wasserstein2_edge_edge_perp(const Segment &e0, const Segment &e1)
{
  std::cout << "The Wasserstein2 distance between segements " << pretty_print(e0) << " and " << pretty_print(e1) << " is" << std::endl;

  const double c = sqrt(1/3.);

  // vectors between some pairs of the four points
  // ensure we match the closest ones
  auto ds = e0.source() - e1.source();
  auto dt = e0.target() - e1.target();
  // for the perp case, the dot product of these is zero
  auto d0 = e0.source() - e0.target();
  auto d1 = e1.source() - e1.target();
  
  // debug
  std::cout << "ds: " << ds << " length:" << sqrt( CGAL::to_double( ds.squared_length() ) ) << std::endl;;
  std::cout << "dt: " << dt << " length:" << sqrt( CGAL::to_double( dt.squared_length() ) ) << std::endl;
  std::cout << "d0: " << d0 << " length:" << sqrt( CGAL::to_double( d0.squared_length() ) ) << std::endl;
  std::cout << "d1: " << d1 << " length:" << sqrt( CGAL::to_double( d1.squared_length() ) ) << std::endl;
  
  auto three_dots = ds * ds + dt * dt + d0 * d1;
  const double W2_analytic = c * sqrt( CGAL::to_double( three_dots ) );
  std::cout << "W2_analytic : " << W2_analytic << std::endl;

  // verify by doing three dot products manually
  if (1)
  {
    K_real ddss(0.), ddtt(0), dd01(0);
    for (int i = 0; i < ds.dimension(); ++i)
    {
      ddss  += ds.cartesian(i) * ds.cartesian(i);
      ddtt  += dt.cartesian(i) * dt.cartesian(i);
      dd01 +=  d0.cartesian(i) * d1.cartesian(i);
    }
    std::cout << "dss: " << CGAL::to_double(ddss) << " ds*ds " << CGAL::to_double(ds*ds) << std::endl;
    
    const double W2_analytic_manual = c * sqrt( CGAL::to_double( ddss + ddtt + dd01 ) );
    std::cout << "W2_manual   : " << W2_analytic_manual << std::endl;
  }

  
  // verification against numeric integration
  if (1)
  {
    // numeric
    //    for (int n = 10; n < 123456; n *= 10)
    for (int n = 1000; n < 1234; n *= 10)
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

double wasserstein2_edge_edge(const Segment &e0, const Segment &e1)
{
  // std::cout << "Computing the Wasserstein distance between segements " << pretty_print(e0) << " and " << pretty_print(e1) << std::endl;

  
  // dumb brute force way
  // each integration point of one segment is mapped to an integration point in another

  // n = number of integration points,
  // if using exact,   n <= 7 for fast, n = 8 involves 15 second wait
  // if using Cartesian, n=8 is instant, n = 10 takes 30 seconds,
  const int n=6;
  assert( n > 1 );

  // data for permutations
  // numbers from 0 to n-1
  std::vector<int> ind(n);
  std::iota(ind.begin(), ind.end(), 0);
  std::vector<int> best_ind(ind);
  std::vector<int> worst_ind(ind);

  // consider all possible pairings of integration points
  // go through all the points of e0 in order
  // go through all the points of e1 in all possible permutations

  const double min_d2d = std::numeric_limits<double>::max();
  K_real min_d2( min_d2d );
  K_real worst_d2( 0. );

  // for all permutations
  do
  {
    
    K_real d2(0.);
    
    // all points of e0
    for (auto i = 0; i < n; ++i)
    {
      Point p0 = frac( e0, i, n );
      Point p1 = frac( e1, ind[i], n );
    
      d2 += CGAL::squared_distance( p0, p1 );  // W2
      // d2 += sqrt( CGAL::to_double( CGAL::squared_distance( p0, p1 ) ) ); // W1

//      std::cout << "i: " << i << " p:" << p0 << " to j:" << ind[i] << " q:" << p1 << " d^2:" << d2 << " == " << d2_float << std::endl;
    }
//    std::cout << "d2:" << d2 << std::endl;
//     std::cout << "  " << d2;
    
    // d2 = sum of distances^2 for this permutation
    // best so far?
    if ( d2 < min_d2 )
    {
      best_ind = ind;
      min_d2 = d2;
//      std::cout << "new min :" << min_d2 << std::endl;
    }
    if ( d2 > worst_d2 )
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
  if (0)
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
  if (0)
  {
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,1), Point(0,3)) );      // 2.43 W2
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,0), Point(0,2)) );      // 1.71
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-0.3), Point(0,1.7)) ); // 1.55
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-0.5), Point(0,1.5)) ); // 1.47
    wasserstein2_edge_edge_perp( Segment(p0,p1), Segment( Point(0,-1), Point(0,1)) );     // 1.39
  }
  
  // skew
  if (0)
  {
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(-1,1), Point(1,1)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(-1,1), Point(1,2)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(-1,1), Point(1,3)) );

    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,1)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,2)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,1), Point(1.8,3)) );

    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.8,0.3)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.8,0.9)) );
    wasserstein2_edge_edge( Segment(p0,p1), Segment( Point(0.2,-0.3), Point(1.8,1.5)) );
  }
}

int main(int argc, char **argv) {
  
  // Test Wasserstein distance between two segments
  if (0)
  {
    wasserstein_edge_edge_test();
    return 0;
  }
  
  // Test hot energy between triangle and triangle* == point
  DT dt;
  const int num_points = 20;
  generate_rand_dt(num_points, dt);
  std::cout << "Generated the Delaunay Triangulation of "
            << num_points << " points" << std::endl;
  std::cout << "Resulting in " << dt.number_of_faces()
            << " faces" << std::endl;
  int face_idx = 1;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end();
      face_itr++, face_idx++) {
    Face face = *face_itr;
    std::cout << "Face " << face_idx << " : ";
    for(int i = 0; i < tri_verts; i++) {
      std::cout << " ( " << face.vertex(i)->point()
                << " )  ";
    }
    std::cout << std::endl;
  }
  K::RT energy = hot_energy(dt);
  std::cout << "Mesh energy: " << energy << std::endl;
  return 0;
}
