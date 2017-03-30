// hot.cpp

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>

#include <polynomial.hpp>

#include <iostream>

constexpr const int dims = 2;

using K = CGAL::Cartesian<double>;
// real type underlying K
using K_real = K::RT;

using DT = CGAL::Delaunay_triangulation_2<K>;
using Face = DT::Face;
using Point = DT::Point;
using Vertex = DT::Vertex;

using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

constexpr const int tri_verts = 3;
constexpr const int tri_edges = 3;

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
 * piecewise integral of the Wasserstein distance with
 * k */
boost::variant<std::array<Triangle, 2>,
               std::array<Triangle, 1> >
integral_bounds(const Triangle &face) {
  std::array<Point, tri_verts> verts = {
      face.vertex(0), face.vertex(1), face.vertex(2)};
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

Triangle face_to_tri(const Face &face) {
  return Triangle(face.vertex(0)->point(),
                  face.vertex(1)->point(),
                  face.vertex(2)->point());
}

std::string pretty_print(const Segment &e) {
  std::stringstream ss;
  ss << "(" << e.source() << ")-(" << e.target() << ")";
  return ss.str();
}

Point frac(const Segment &e, int i, int n) {
  return e.source() + (i / (double)(n - 1)) * e.to_vector();
}

double wasserstein2_edge_edge(const Segment &e0,
                              const Segment &e1) {
  // std::cout << "Computing the Wasserstein distance
  // between segements " << pretty_print(e0) << " and " <<
  // pretty_print(e1) << std::endl;

  // dumb brute force way
  // each integration point of one segment is mapped to an
  // integration point in another

  // n = number of integration points,
  // if using exact,   n <= 7 for fast, n = 8 involves 15
  // second wait
  // if using Cartesian, n=8 is instant, n = 10 takes 30
  // seconds,
  const int n = 6;
  assert(n > 1);

  // data for permutations
  // numbers from 0 to n-1
  std::vector<int> ind(n);
  std::iota(ind.begin(), ind.end(), 0);
  std::vector<int> best_ind(ind);
  std::vector<int> worst_ind(ind);

  // consider all possible pairings of integration points
  // go through all the points of e0 in order
  // go through all the points of e1 in all possible
  // permutations

  const double min_d2d = std::numeric_limits<double>::max();
  K_real min_d2(min_d2d);
  K_real worst_d2(0.);

  // for all permutations
  do {
    K_real d2(0.);

    // all points of e0
    for(auto i = 0; i < n; ++i) {
      Point p0 = frac(e0, i, n);
      Point p1 = frac(e1, ind[i], n);

      d2 += CGAL::squared_distance(p0, p1);  // W2
      // d2 += sqrt( CGAL::to_double(
      // CGAL::squared_distance( p0, p1 ) ) ); // W1

      //      std::cout << "i: " << i << " p:" << p0 << " to
      //      j:" << ind[i] << " q:" << p1 << " d^2:" << d2
      //      << " == " << d2_float << std::endl;
    }
    //    std::cout << "d2:" << d2 << std::endl;
    //     std::cout << "  " << d2;

    // d2 = sum of distances^2 for this permutation
    // best so far?
    if(d2 < min_d2) {
      best_ind = ind;
      min_d2 = d2;
      //      std::cout << "new min :" << min_d2 <<
      //      std::endl;
    }
    if(d2 > worst_d2) {
      worst_ind = ind;
      worst_d2 = d2;
      //      std::cout << "new worst :" << worst_d2 <<
      //      std::endl;
    }

  } while(std::next_permutation(ind.begin(), ind.end()));
  // if permutation speed becomes important, try
  // http://howardhinnant.github.io/combinations.html

  auto min_d = sqrt(CGAL::to_double(min_d2 / n));      // W2
  auto worst_d = sqrt(CGAL::to_double(worst_d2 / n));  // W2
  // auto min_d = CGAL::to_double(min_d2) / (double) n; //
  // W1
  // auto worst_d = CGAL::to_double(worst_d2) / (double) n;
  // // W1

  std::cout
      << "\nThe Wasserstein2 distance between segements "
      << pretty_print(e0) << " and " << pretty_print(e1)
      << " was " << min_d << std::endl;
  std::cout << "From permutation  :";
  for(auto p : best_ind) std::cout << " " << p;
  std::cout << "\nWorst permutation :";
  for(auto p : worst_ind) std::cout << " " << p;
  std::cout << " gave distance " << worst_d << ". ";

  std::cout << std::endl << std::endl;

  return min_d;
}
