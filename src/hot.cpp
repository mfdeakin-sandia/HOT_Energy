
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>

#include <cmath>
#include <random>

#include <iostream>

using RNG = std::mt19937_64;

using Real = float;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using DT = CGAL::Delaunay_triangulation_2<K>;
using Face = DT::Face;
using Point = DT::Point;

constexpr const int dims = 2;

constexpr const int tri_verts = 3;

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
  Real x_c[dims] = {0.0, 0.0};
  for(int i = 0; i < tri_verts; i++) {
    for(int j = 0; j < dims; j++) {
      x_c[j] +=
          CGAL::to_double(face.vertex(i)->point()[j] / 3.0);
    }
  }
  return Point(x_c[0], x_c[1]);
}

/* Computes the k Wasserstein distance of the triangular
 * face to it's centroid */
template <int k>
Real triangle_w(const Face &face);

template <>
Real triangle_w<1>(const Face &face) {
  Point centroid = triangle_centroid(face);
  // Set to NaN until implemented
  Real distance = 0.0 / 0.0;
  return distance;
}

/* See "HOT: Hodge-Optimized Triangulations" for details on
 * the energy functional */
Real hot_energy(const DT &dt) {
  Real energy = 0.0;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end(); face_itr++) {
    energy +=
        triangle_w<1>(*face_itr) * triangle_area(*face_itr);
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
