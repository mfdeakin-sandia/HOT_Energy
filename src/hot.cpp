// hot.cpp

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <iostream>

#include <polynomial.hpp>

using RNG = std::mt19937_64;

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

constexpr const float min_pos = 10.0;
constexpr const float max_pos = 20.0;

void generate_rand_dt(int num_points, DT &dt) {
  std::random_device rd;
  RNG engine(rd());
  std::uniform_real_distribution<double> genPos(min_pos,
                                                max_pos);
  for(int i = 0; i < num_points; i++) {
    Point pt(genPos(engine), genPos(engine));
    dt.insert(pt);
  }
}

/* Computes the centroid of the triangle */
Point triangle_centroid(const Triangle &face) {
  return CGAL::centroid(face.vertex(0), face.vertex(1),
                        face.vertex(2));
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

template <typename T, int k>
class triangle_w_helper;

template <typename T>
class triangle_w_helper<T, 2> {
 public:
  K_real operator()(T &bounds, const Point &centroid) {
    K_real integral = 0;
    for(auto area : bounds) {
      /* Start by computing the lines used for boundaries
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
      K_real width = verts[w_idx][0] -
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
      K_real slope_1 = dir_1[1] / dir_1[0];
      K_real slope_2 = dir_2[1] / dir_2[0];
      if((w_idx == 0 && slope_1 < slope_2) ||
         (w_idx == 2 && slope_1 > slope_2)) {
        std::swap(slope_1, slope_2);
      }
      Numerical::Polynomial<K_real, 1, 1> bound_1;
      bound_1.coeff(1) = slope_1;
      bound_1.coeff(0) =
          -verts[w_idx][0] * slope_1 + verts[w_idx][1];
      Numerical::Polynomial<K_real, 1, 1> bound_2;
      bound_2.coeff(1) = slope_2;
      bound_2.coeff(0) =
          -verts[w_idx][0] * slope_2 + verts[w_idx][1];

      Numerical::Polynomial<K_real, 1, 2> initial_x_root(
          (Tags::Zero_Tag()));
      initial_x_root.coeff(1, 0) = 1;
      initial_x_root.coeff(0, 0) = -centroid[0];

      Numerical::Polynomial<K_real, 1, 2> initial_y_root(
          (Tags::Zero_Tag()));
      initial_y_root.coeff(0, 1) = 1;
      initial_y_root.coeff(0, 0) = -centroid[1];

      Numerical::Polynomial<K_real, 2, 2> initial =
          initial_x_root * initial_x_root +
          initial_y_root * initial_y_root;

      auto y_int = initial.integrate(1);
      Numerical::Polynomial<K_real, 3, 1> upper =
          y_int.var_sub(1, bound_1);
      // 3.55556 y + -1.33333 y ^ 2 + 0.333333 y ^ 3 +
      // -2.66667 x y + x ^ 2 y
      // y = y_1 = -x + 3
      // -1.33333 + -3.55556 x + 4.33333 x ^ 2 + -1 x ^ 3
      Numerical::Polynomial<double, 3, 1> lower =
          y_int.var_sub(1, bound_2);
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
K_real triangle_w(const Triangle &face);

template <>
K_real triangle_w<2>(const Triangle &tri) {
  Point centroid = triangle_centroid(tri);
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  // Set to NaN until implemented
  K_real distance = 0.0 / 0.0;
  if(bounds.which() == 0) {
    // std::array<Triangle, 2>
    auto area =
        boost::get<std::array<Triangle, 2> >(bounds);
    return triangle_w_helper<std::array<Triangle, 2>, 2>()(
        area, centroid);
  } else {
    // std::array<Triangle, 1>
    auto area =
        boost::get<std::array<Triangle, 1> >(bounds);
    return triangle_w_helper<std::array<Triangle, 1>, 2>()(
        area, centroid);
  }
}

Triangle face_to_tri(const Face &face) {
  return Triangle(face.vertex(0)->point(),
                  face.vertex(1)->point(),
                  face.vertex(2)->point());
}

/* See "HOT: Hodge-Optimized Triangulations" for details on
 * the energy functional */
K_real hot_energy(const DT &dt) {
  K_real energy = 0;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end(); face_itr++) {
    Triangle tri = face_to_tri(*face_itr);
    K_real area = tri.area();
    K_real wasserstein = triangle_w<2>(tri);
    energy += wasserstein * area;
  }
  return energy;
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

void wasserstein_edge_edge_test() {
  // create two test edges
  Point p0(0, 0);
  Point p1(2, 0);
  Point p2(1, 0);
  Point p3(1, 2);

  // use P10--P11 for both, so distance should be zero
  Segment e0(p0, p1);
  Segment e1(p0, p1);

  // compute distance
  //  wasserstein2_edge_edge( Segment(p0,p1), Segment(p0,p1)
  //  );
  //  wasserstein2_edge_edge( Segment(p0,p1), Segment(p1,p0)
  //  );

  // perp edges at midpoint
  // of uniform length
  // appears that W2 is invariant to permutation of matching
  // up integration points!
  // W1 order matters, best is a large spread, match close
  // points together and match far points together
  // W1 best is interlaced, alternating points from each
  // side, *not* 0-5 matched with 0-5
  if(0) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, 1), Point(1, 3)));  // 2.22  W2
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, 0), Point(1, 2)));  // 1.39
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -0.3), Point(1, 1.7)));  // 1.19
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -0.5), Point(1, 1.5)));  // 1.09
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -1), Point(1, 1)));  // 0.966
  }

  // perp edges, not at midpoint
  // appears that W2 is invariant to permutation of matching
  // up integration points!
  // W1 is as in the above case
  if(0) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, 1), Point(0, 3)));  // 2.43 W2
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, 0), Point(0, 2)));  // 1.71
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -0.3), Point(0, 1.7)));  // 1.55
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -0.5), Point(0, 1.5)));  // 1.47
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -1), Point(0, 1)));  // 1.39
  }

  // skew
  if(1) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 1)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 2)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 3)));

    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 1)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 2)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 3)));

    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 0.3)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 0.9)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 1.5)));
  }
}

void test_tri_w2() {
  Triangle tri(Point(2, 1), Point(3, 1), Point(2, 2));
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  assert(bounds.which() == 1);
  auto area = boost::get<std::array<Triangle, 1> >(bounds);

  std::cout
      << std::endl
      << "Computing Wasserstein distance to the dual point"
      << std::endl
      << std::endl;
  K_real value = triangle_w<2>(tri);
  std::cout << "Unit Right Triangle expected Wasserstein "
               "Distance: 5/90 (~0.0555555). Calculated "
               "Distance: "
            << value << std::endl
            << std::endl;
}

int main(int argc, char **argv) {
  // Test Wasserstein distance between two segments
  if(0) {
    wasserstein_edge_edge_test();
    test_tri_w2();
    return 0;
  }

  // Test hot energy between triangle and triangle* == point

  DT dt;
  const int num_points = 30;
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
  K_real energy = hot_energy(dt);
  std::cout << "Mesh energy: " << energy << std::endl;
  return 0;
}
