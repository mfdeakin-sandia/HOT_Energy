
#ifndef _HOT_OPTIMIZED_MESH_HPP_
#define _HOT_OPTIMIZED_MESH_HPP_

#include <OptimizedMesh.hpp>

#include <CGAL/Delaunay_Triangulation_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>

#include <boost/variant/variant.hpp>

#include <polynomial.hpp>

#include <array>

using K = CGAL::Cartesian<real>;
using DT = CGAL::Delaunay_triangulation_2<K>;

using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

template <typename real> class HotOptimizedMesh : OptimizedMesh<DT, real> {
public:
private:
  virtual real energy() {
    real energy = 0.0;
    for (auto face_itr = dt.finite_faces_begin();
         face_itr != dt.finite_faces_end(); face_itr++) {
      Triangle tri = face_to_tri(*face_itr);
      energy += triangle_energy(tri);
    }
    return energy;
  }

  virtual real triangle_energy(const Triangle &tri) {
    Point centroid = triangle_centroid(tri);
    boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1> > bounds =
        integral_bounds(tri);
    real energy = 0.0;
    if (bounds.which() == 0) {
      // std::array<Triangle, 2>
      auto area = boost::get<std::array<Triangle, 2> >(bounds);
      for (Triangle tri : area) {
        energy += triangle_wasserstein_2(tri, centroid);
      }
    } else {
      // std::array<Triangle, 1>
      auto area = boost::get<std::array<Triangle, 1> >(bounds);
      for (Triangle tri : area) {
        energy += triangle_wasserstein_2(tri, centroid);
      }
    }
    return energy;
  }

  real triangle_wasserstein_2(const Triangle &area, const Point &centroid)
  {
    /* Start by computing the lines used for boundaries
     * and the width.
     */
    std::array<Point, tri_verts> verts = { area.vertex(0), area.vertex(1),
                                           area.vertex(2) };
    order_points(verts);

    /* Find the index of the non-vertical point
     * It's guaranteed to be either the leftmost or
     * rightmost point
     */
    const int w_idx = (verts[0][0] != verts[1][0]) ? 0 : 2;
    // Compute the signed width
    K_real width = verts[w_idx][0] - verts[(w_idx + 1) % tri_verts][0];
    /* Compute the bounding lines
     * x = ax t + bx -> t = (x - bx) / ax
     * y = ay t + by = (ay / ax) x + (-ay bx / ax + by)
     *
     * bx = verts[w_idx][0], by = verts[w_idx][1]
     */
    Vector dir_1 =
        Line(verts[w_idx], verts[(w_idx + 1) % tri_verts]).to_vector();
    Vector dir_2 =
        Line(verts[w_idx], verts[(w_idx + 2) % tri_verts]).to_vector();
    K_real slope_1 = dir_1[1] / dir_1[0];
    K_real slope_2 = dir_2[1] / dir_2[0];
    if ((w_idx == 0 && slope_1 < slope_2) ||
        (w_idx == 2 && slope_1 > slope_2)) {
      std::swap(slope_1, slope_2);
    }
    Numerical::Polynomial<K_real, 1, 1> bound_1;
    bound_1.coeff(1) = slope_1;
    bound_1.coeff(0) = -verts[w_idx][0] * slope_1 + verts[w_idx][1];
    Numerical::Polynomial<K_real, 1, 1> bound_2;
    bound_2.coeff(1) = slope_2;
    bound_2.coeff(0) = -verts[w_idx][0] * slope_2 + verts[w_idx][1];

    Numerical::Polynomial<K_real, 1, 2> initial_x_root((Tags::Zero_Tag()));
    initial_x_root.coeff(1, 0) = 1;
    initial_x_root.coeff(0, 0) = -centroid[0];

    Numerical::Polynomial<K_real, 1, 2> initial_y_root((Tags::Zero_Tag()));
    initial_y_root.coeff(0, 1) = 1;
    initial_y_root.coeff(0, 0) = -centroid[1];

    Numerical::Polynomial<K_real, 2, 2> initial =
        initial_x_root * initial_x_root + initial_y_root * initial_y_root;

    auto y_int = initial.integrate(1);
    Numerical::Polynomial<K_real, 3, 1> upper = y_int.var_sub(1, bound_1);
    Numerical::Polynomial<double, 3, 1> lower = y_int.var_sub(1, bound_2);
    auto y_bounded = upper + -lower;

    auto x_int = y_bounded.integrate(0);

    auto left = x_int.slice(0, verts[w_idx][0]).coeff(0);
    auto right = x_int.slice(0, verts[(w_idx + 1) % tri_verts][0]).coeff(0);
    real integral = std::abs(left + -right);
    return integral;
  }
};

#endif
