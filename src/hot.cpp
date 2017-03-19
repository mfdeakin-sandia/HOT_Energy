
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <iostream>

#include <polynomial.hpp>

using RNG = std::mt19937_64;

using Real = float;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using DT = CGAL::Delaunay_triangulation_2<K>;
using Face = DT::Face;
using Point = DT::Point;

using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

constexpr const int dims = 2;

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

int main(int argc, char **argv) {
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
