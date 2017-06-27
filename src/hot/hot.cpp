// hot.cpp

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <list>

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
  if (verts[0][0] > verts[1][0]) {
    std::swap(verts[0], verts[1]);
  }
  if (verts[1][0] > verts[2][0]) {
    std::swap(verts[2], verts[1]);
    if (verts[0][0] > verts[1][0]) {
      std::swap(verts[0], verts[1]);
    }
  }
}

/* Computes 1 or 2 triangles which can be used for the
 * piecewise integral of the Wasserstein distance with
 * k */
boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1> >
integral_bounds(const Triangle &face) {
  std::array<Point, tri_verts> verts = { face.vertex(0), face.vertex(1),
                                         face.vertex(2) };
  order_points(verts);
  if (verts[0][0] == verts[1][0] || verts[1][0] == verts[2][0]) {
    // Return a single triangle in this case
    std::array<Triangle, 1> bounds;
    bounds[0] = Triangle(verts[0], verts[1], verts[2]);
    return boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1> >(
        bounds);
  } else {
    const Line vertical(verts[1], verts[1] + Vector(0, 1));
    const Line base(verts[0], verts[2]);
    // This intersection will exist for all non-degenerate
    // triangles
    auto int_vert = CGAL::intersection(vertical, base);
    assert(int_vert.is_initialized());

    std::array<Triangle, 2> bounds;
    bounds[0] = Triangle(verts[0], verts[1], boost::get<Point>(int_vert.get()));
    bounds[1] = Triangle(verts[2], verts[1], boost::get<Point>(int_vert.get()));
    return boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1> >(
        bounds);
  }
}

Triangle face_to_tri(const Face &face) {
  return Triangle(face.vertex(0)->point(), face.vertex(1)->point(),
                  face.vertex(2)->point());
}

Point triangle_circumcenter(const Triangle &face) {
  return CGAL::circumcenter(face.vertex(0), face.vertex(1), face.vertex(2));
}

std::list<DT::Vertex_handle> internal_vertices(const DT &dt) {
  std::list<DT::Vertex_handle> verts;
  std::set<DT::Vertex_handle> seen;
  // First mark the external vertices as having been seen
  for (auto vert_itr = dt.incident_vertices(dt.infinite_vertex());
       seen.count(vert_itr) == 0; vert_itr++) {
    seen.insert(vert_itr);
  }
  // Now iterator over all of the vertices,
  // storing them if we haven't seen them before
  for (auto vert_itr = dt.finite_vertices_begin();
       vert_itr != dt.finite_vertices_end(); vert_itr++) {
    if (seen.count(vert_itr) == 0) {
      seen.insert(vert_itr);
      verts.push_back(vert_itr);
    }
  }
  return verts;
}
