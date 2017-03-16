
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <iostream>

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

Real triangle_area(const Face &face) {
  Point origin = face.vertex(0)->point();
  Point sides[tri_verts - 1];
  for(int i = 0; i < tri_verts - 1; i++) {
    sides[i] =
        Point(face.vertex(i)->point()[0] - origin[0],
              face.vertex(i)->point()[1] - origin[1]);
  }
  Real area = std::fabs(CGAL::to_double(
                  sides[0][0] * sides[1][1] -
                  sides[1][0] * sides[0][1])) /
              2.0;
  return area;
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

/* Computes 1 or 2 triangles which can be used for the
 * piecewise integral of the Wasserstein distance with even
 * k */
boost::variant<std::array<Triangle, 2>,
               std::array<Triangle, 1> >
even_integral_bounds(const Face &face) {
  Point verts[tri_verts] = {face.vertex(0)->point(),
                            face.vertex(1)->point(),
                            face.vertex(2)->point()};
  if(verts[0][0] > verts[1][0]) {
    std::swap(verts[0], verts[1]);
  }
  if(verts[1][0] > verts[2][0]) {
    std::swap(verts[2], verts[1]);
  }
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
  Real operator()(T &bounds) {
    Real integral = 0.0 / 0.0;
    for(auto area : bounds) {
    }
    return integral;
  }
};

/* Computes the k Wasserstein distance of the triangular
 * face to it's centroid */
template <int k>
Real triangle_w(const Face &face);

template <>
Real triangle_w<2>(const Face &face) {
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
Real hot_energy(const DT &dt) {
  Real energy = 0.0;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end(); face_itr++) {
    energy +=
        triangle_w<2>(*face_itr) * triangle_area(*face_itr);
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
  std::cout << "Mesh energy: " << hot_energy(dt)
            << std::endl;
  return 0;
}
