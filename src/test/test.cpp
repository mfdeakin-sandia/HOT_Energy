
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>

#include <iomanip>
#include <iostream>

#include "array.hpp"
#include "hot.hpp"
#include "ply_writer.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using RNG = std::mt19937_64;

template <typename A>
void permute_helper(A &avail_items, int choice_index, int &permutation,
                    std::function<void(const A &, int)> &f) {
  if (choice_index > 0) {
    permute_helper(avail_items, choice_index - 1, permutation, f);
    for (int i = 0; i < choice_index; i++) {
      std::swap(avail_items[i], avail_items[choice_index]);
      permute_helper(avail_items, choice_index - 1, permutation, f);
    }
  } else {
    f(avail_items, permutation);
    permutation++;
  }
}

template <typename A>
void permute(A items, std::function<void(const A &, int)> &f) {
  int permutation = 0;
  permute_helper(items, A::size() - 1, permutation, f);
}

TEST_CASE("Unit Right Triangle", "[HOT]") {
  constexpr const int num_permutations = CTMath::partialFactorial(1, tri_verts);
  constexpr const int num_shifts = 10;
  constexpr const int num_tris = num_shifts * num_permutations;
  // TODO: Reduce the large relative error of the
  // Wasserstein calculation
  // 20 bits of double precision's 53 is too much!
  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 1048376.0;

  std::random_device rd;
  RNG engine(rd());
  constexpr const double min_shift = -10.0;
  constexpr const double max_shift = 9.0;
  // Use float to reduce the number of bits required to
  // represent the final number
  std::uniform_real_distribution<float> genShift(min_shift, max_shift);

  Array<Point, 3> points = {Point(0, 0), Point(1, 0), Point(0, 1)};
  Array<Triangle, num_tris> tri;
  std::function<void(const Array<Point, 3> &, int)> filler(
      [&](const Array<Point, 3> &pt_permute, int permutation) {
        for (int i = 0; i < num_shifts; i++) {
          double dx = genShift(engine);
          double dy = genShift(engine);
          Array<Point, 3> test_point = {
              Point(pt_permute[0][0] + dx, pt_permute[0][1] + dy),
              Point(pt_permute[1][0] + dx, pt_permute[1][1] + dy),
              Point(pt_permute[2][0] + dx, pt_permute[2][1] + dy)};
          tri[i * num_permutations + permutation] =
              Triangle(test_point[0], test_point[1], test_point[2]);
        }
      });
  permute(points, filler);

  SECTION("Split") {
    for (int i = 0; i < num_tris; i++) {
      boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1>> bounds =
          integral_bounds(tri[i]);
      REQUIRE(bounds.which() == 1);
    }
  }
  SECTION("Wasserstein 2") {
    // The correct answer is (5 / 60)
    for (int i = 0; i < num_tris; i++) {
      K_real energy = triangle_w<2>(tri[i]);
      REQUIRE(std::abs(energy * 12 - 1) < max_rel_error);
    }
  }
}

TEST_CASE("Unit Equilateral Triangle", "[HOT]") {
  constexpr const int num_permutations = CTMath::partialFactorial(1, tri_verts);
  constexpr const int num_shifts = 10;
  constexpr const int num_tris = num_shifts * num_permutations;
  // TODO: Reduce the large relative error of the
  // Wasserstein calculation
  // 21 bits of double precision's 53 is too much!
  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 2097152.0;

  std::random_device rd;
  RNG engine(rd());
  constexpr const double min_shift = -10.0;
  constexpr const double max_shift = 9.0;
  // Use float to reduce the number of bits required to
  // represent the final number
  std::uniform_real_distribution<float> genShift(min_shift, max_shift);

  Array<Point, 3> points = {Point(0, 0), Point(1, 0),
                            Point(0.5, std::sqrt(3) / 2.0)};
  Array<Triangle, num_tris> tri;
  std::function<void(const Array<Point, 3> &, int)> filler(
      [&](const Array<Point, 3> &pt_permute, int permutation) {
        for (int i = 0; i < num_shifts; i++) {
          double dx = genShift(engine);
          double dy = genShift(engine);
          Array<Point, 3> test_point = {
              Point(pt_permute[0][0] + dx, pt_permute[0][1] + dy),
              Point(pt_permute[1][0] + dx, pt_permute[1][1] + dy),
              Point(pt_permute[2][0] + dx, pt_permute[2][1] + dy)};
          tri[i * num_permutations + permutation] =
              Triangle(test_point[0], test_point[1], test_point[2]);
        }
      });
  permute(points, filler);

  SECTION("Split") {
    for (int i = 0; i < num_tris; i++) {
      boost::variant<std::array<Triangle, 2>, std::array<Triangle, 1>> bounds =
          integral_bounds(tri[i]);
      REQUIRE(bounds.which() == 0);
    }
  }
  SECTION("Wasserstein 2") {
    for (int i = 0; i < num_tris; i++) {
      K_real energy = triangle_w<2>(tri[i]);
      REQUIRE(std::abs(energy - 0.036084391824351578232) <
              0.036084391824351578232 * max_rel_error);
    }
  }
}

TEST_CASE("Single Point Mesh Gradient Descent Local Minimum", "[HOT]") {
  using K = CGAL::Cartesian<double>;
  using K_real = K::RT;
  using DT = CGAL::Delaunay_triangulation_2<K>;

  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 32.0;

  constexpr const int num_bounds = 4;
  const double x_bounds[] = {-std::sqrt(0.5), -std::sqrt(0.5), std::sqrt(0.5),
                             std::sqrt(0.5)};
  const double y_bounds[] = {-std::sqrt(0.5), std::sqrt(0.5), std::sqrt(0.5),
                             -std::sqrt(0.5)};
  constexpr const double initial_internal_x = 0.0;
  constexpr const double initial_internal_y = 0.0;

  DT dt;
  for (int i = 0; i < num_bounds; i++) {
    DT::Point pt(x_bounds[i], y_bounds[i]);
    dt.insert(pt);
  }

  dt.insert(DT::Point(initial_internal_x, initial_internal_y));

  std::list<DT::Vertex_handle> internal_verts = internal_vertices(dt);

  SECTION("Internal Vertices") {
    // Verify there's only 1 internal vertex
    REQUIRE(internal_verts.size() == 1);
    // Vertify that the vertex is the one we expect
    DT::Vertex_handle vertex = internal_verts.front();
    REQUIRE(vertex->point()[0] == initial_internal_x);
    REQUIRE(vertex->point()[1] == initial_internal_y);
  }

  SECTION("Vertex Gradient") {
    // This is a local minimum, so the energy should be iso
    std::vector<finite_diffs> f_diffs = compute_gradient<2>(dt, internal_verts);
    REQUIRE(f_diffs.size() == 1);
    // This is the sum of 4 right triangles of unit side length
    // The energy is then 4 * 1/2 * 5/60 = 1/6
    REQUIRE(std::abs(f_diffs[0].center * 6.0 - 1.0) <= max_rel_error);
    // Because we're at a local minimum, the figure is symmetric,
    // and we're moving the same in either direction,
    // dx_plus == dx_minus == dy_plus == dy_minus
    REQUIRE(std::abs(f_diffs[0].dx_plus - f_diffs[0].dx_minus) <=
            max_rel_error);
    REQUIRE(std::abs(f_diffs[0].dy_plus - f_diffs[0].dy_minus) <=
            max_rel_error);
  }

  SECTION("Vertex Gradient Descent") {
    DT optimized = hot_optimize<2>(dt);
    std::list<DT::Vertex_handle> optimized_verts = internal_vertices(optimized);
    REQUIRE(optimized_verts.size() == 1);
    DT::Vertex_handle vertex = optimized_verts.front();
    REQUIRE(std::abs(vertex->point()[0] - initial_internal_x) <= max_rel_error);
    REQUIRE(std::abs(vertex->point()[1] - initial_internal_y) <= max_rel_error);
  }
}

TEST_CASE("Single Point Mesh Gradient Descent", "[HOT]") {
  using K = CGAL::Cartesian<double>;
  using K_real = K::RT;
  using DT = CGAL::Delaunay_triangulation_2<K>;

  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 32.0;

  constexpr const int num_bounds = 3;
  const double x_bounds[] = {-1.0, -1.0, 1.0};
  const double y_bounds[] = {-1.0, 1.0, 0.0};
  constexpr const double initial_internal_x = 0.875;
  constexpr const double initial_internal_y = 0.0;

  DT dt;
  for (int i = 0; i < num_bounds; i++) {
    DT::Point pt(x_bounds[i], y_bounds[i]);
    dt.insert(pt);
  }

  dt.insert(DT::Point(initial_internal_x, initial_internal_y));

  std::list<DT::Vertex_handle> internal_verts = internal_vertices(dt);

  SECTION("Internal Vertices") {
    // Verify there's only 1 internal vertex
    REQUIRE(internal_verts.size() == 1);
    // Vertify that the vertex is the one we expect
    DT::Vertex_handle vertex = internal_verts.front();
    REQUIRE(vertex->point()[0] == initial_internal_x);
    REQUIRE(vertex->point()[1] == initial_internal_y);
  }

  SECTION("Vertex Gradient") {
    // This is a local minimum, so the energy should be iso
    std::vector<finite_diffs> f_diffs = compute_gradient<2>(dt, internal_verts);
    REQUIRE(f_diffs.size() == 1);

    // This figure was constructed only so that the gradient will be in the
    // negative x direction, and 0 in the y direction
    REQUIRE(f_diffs[0].dx_minus < f_diffs[0].dx_plus - max_rel_error);
    REQUIRE(std::abs(f_diffs[0].dy_plus - f_diffs[0].dy_minus) <=
            max_rel_error);
  }

  SECTION("Vertex Gradient Descent") {
    DT optimized = hot_optimize<2>(dt);
    std::list<DT::Vertex_handle> optimized_verts = internal_vertices(optimized);
    REQUIRE(optimized_verts.size() == 1);
    DT::Vertex_handle vertex = optimized_verts.front();
    REQUIRE(vertex->point()[0] < initial_internal_x);
    REQUIRE(std::abs(vertex->point()[1] - initial_internal_y) <= max_rel_error);
  }
}

TEST_CASE("Two Point Mesh Gradient Descent", "[HOT]") {
  using K = CGAL::Cartesian<double>;
  using K_real = K::RT;
  using DT = CGAL::Delaunay_triangulation_2<K>;

  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 32.0;

  constexpr const int num_bounds = 6;
  const double x_bounds[] = {-1.5, -0.5, 0.5, 1.5, 0.5, -0.5};
  const double y_bounds[] = {-0.5, 0.5, 1.5, 0.5, -0.5, -1.5};

  constexpr const int num_internal = 2;
  constexpr const double initial_internal_x[] = {-0.5, 0.5};
  constexpr const double initial_internal_y[] = {-0.5, 0.5};

  DT dt;
  for (int i = 0; i < num_bounds; i++) {
    DT::Point pt(x_bounds[i], y_bounds[i]);
    dt.insert(pt);
  }

  for (int i = 0; i < num_internal; i++) {
    dt.insert(DT::Point(initial_internal_x[i], initial_internal_y[i]));
  }

  std::list<DT::Vertex_handle> internal_verts = internal_vertices(dt);

  SECTION("Internal Vertices") {
    // Verify there's only 1 internal vertex
    REQUIRE(internal_verts.size() == 2);
    // Vertify that the vertex is the one we expect
    DT::Vertex_handle vertex = internal_verts.front();
    REQUIRE(vertex->point()[0] == initial_internal_x[0]);
    REQUIRE(vertex->point()[1] == initial_internal_y[0]);

    vertex = internal_verts.back();
    REQUIRE(vertex->point()[0] == initial_internal_x[1]);
    REQUIRE(vertex->point()[1] == initial_internal_y[1]);
  }

  SECTION("Vertex Gradient") {
    // This is a local minimum, so the energy should be iso
    std::vector<finite_diffs> f_diffs = compute_gradient<2>(dt, internal_verts);
    REQUIRE(f_diffs.size() == 2);
    // This is the sum of 4 right triangles of unit side length
    // The energy is then 4 * 1/2 * 5/60 = 1/6
    REQUIRE(std::abs(f_diffs[0].center * 6.0 - 1.0) <= max_rel_error);
    REQUIRE(std::abs(f_diffs[1].center * 6.0 - 1.0) <= max_rel_error);
  }
}
