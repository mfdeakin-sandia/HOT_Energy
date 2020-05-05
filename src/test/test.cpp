
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

TEST_CASE("HOT_fig1", "[HOT]")
{
  
  // read file
  
  // directory path?
  std::string epath;
  auto energy_dir_p = getenv("HOT_Energy");
  // expect something like
  //   /Users/samitch/Documents/repos/PrimalDual/HOT_Energy
  if (energy_dir_p==nullptr)
  {
    std::cout << "Warning, environment variable 'HOT_Energy' should be the file location of your clone of the HOT_Energy github repository, but is undefined" << std::endl;
  }
  else
  {
    epath = energy_dir_p;
  }
  
  // example path
  std::string example_name = "HOT_fig1";
  std::string fpath = epath + "/examples/" + example_name + "/" + example_name;

  // points
  std::vector<Point_2> points;
  {
    std::ifstream p_in(fpath + "_points.txt");
    std::istream_iterator<Point> p_begin(p_in), p_end;
    points.assign(p_begin, p_end);
  }
  
  // weights
  std::vector<double> weights;
  {
    std::ifstream w_in(fpath + "_weights.txt");
    std::istream_iterator<double> w_begin(w_in), w_end;
    weights.assign(w_begin, w_end);
  }
  
  // triangles
  std::vector< int > triangle_points;
  // point 0,1,2 define the first triangle, etc.
  // In the file, indices start at 1, so subtract 1 here
  {
    std::ifstream t_in(fpath + "_triangles.txt");
    std::istream_iterator< int > t_begin(t_in), t_end;
    triangle_points.assign(t_begin, t_end);
  }
  
  for (auto &p : points)
  {
    std::cout << p;
  }
  std::cout << std::endl;
  
  using K = CGAL::Cartesian<double>;
  using K_real = K::RT;

  if (points.empty())
    return;
  
  // find bounding box of points
  double x_bounds[2], y_bounds[2];
  x_bounds[0] = x_bounds[1] = points[0].x();
  y_bounds[0] = y_bounds[1] = points[0].y();
  for (auto &p : points)
  {
    x_bounds[0] = std::min( x_bounds[0], p.x() );
    x_bounds[1] = std::max( x_bounds[1], p.x() );
    y_bounds[0] = std::min( y_bounds[0], p.y() );
    y_bounds[1] = std::max( y_bounds[1], p.y() );
  }

  double x_box[2], y_box[2];
  const double buff = 0.5;
  x_box[0] = x_bounds[0]- buff*(x_bounds[1]-x_bounds[0]);
  x_box[1] = x_bounds[1]+ buff*(x_bounds[1]-x_bounds[0]);
  y_box[0] = y_bounds[0]- buff*(y_bounds[1]-y_bounds[0]);
  y_box[1] = y_bounds[1]+ buff*(y_bounds[1]-y_bounds[0]);

  std::vector< Point_2 > b_points;
  b_points.push_back( DT::Point(x_box[0], y_box[0]) );
  b_points.push_back( DT::Point(x_box[1], y_box[1]) );
  b_points.push_back( DT::Point(x_box[0], y_box[1]) );
  b_points.push_back( DT::Point(x_box[1], y_box[0]) );

  
  // make a Delaunay triangulation of points
  using DT = CGAL::Delaunay_triangulation_2<K>;
  DT dt;

  // insert bounding box points
  // not needed, because it will instead create the point at infinity
  dt.insert( b_points.begin(), b_points.end() );
  
  // insert actual points
  dt.insert(points.begin(),points.end());
  
  // make a weighted DT of points
  using WDT = CGAL::Regular_triangulation_2<K>;
  WDT wdt;
  
  using WP = Regular_triangulation_2::Weighted_point;
  // put points and weights together
  std::vector<WP> wpoints;
  for (size_t i = 0; i < points.size(); ++i)
  {
    wpoints.push_back( WP( points[i], weights[i] ) );
  }
  // insert actual points
  wdt.insert(wpoints.begin(),wpoints.end());
  
  
  // make triangulation of exactly the specified triangles
  // see https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html
  // https://doc.cgal.org/latest/Triangulation_2/index.html#Section_2D_Triangulations_Basic

  // the relevant methods are somewhat "hidden" because these are "advanced" uses which do not guarantee the validity of the triangulation
  // https://doc.cgal.org/latest/TDS_2/classTriangulationDataStructure__2.html
  // see the "create_face" methods
  
  // CGAL::Triangulation_2<K> t; // not this, it does flips and such automagically
  CGAL::Triangulation_data_structure_2<> td;
  // create vertices corresponding to each of the points
  std::vector<CGAL::Triangulation_data_structure_2<>::Vertex_handle> verts;
  for (auto p : points)
  {
    verts.push_back( td.create_vertex() );
  }
  for (auto tpi = triangle_points.begin(); tpi < triangle_points.end(); tpi += 3)
  {
    td.create_face( verts[*tpi], verts[*(tpi+1)], verts[*(tpi+2)] );
    // do we need to tell the faces who their face neighbors are? Or is this sufficient?
  }
  
  // we're on our own for deciding if a triangle is inside out, and calling the corresponding flip functions, etc.
  
}
