
#include <cmath>
#include <limits>

#include <algorithm>

#include <random>

#include <iomanip>
#include <iostream>

#include <array.hpp>
#include <hot.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using RNG = std::mt19937_64;

template <typename A>
void permute_helper(
    A &avail_items, int choice_index, int &permutation,
    std::function<void(const A &, int)> &f) {
  if(choice_index > 0) {
    permute_helper(avail_items, choice_index - 1,
                   permutation, f);
    for(int i = 0; i < choice_index; i++) {
      std::swap(avail_items[i], avail_items[choice_index]);
      permute_helper(avail_items, choice_index - 1,
                     permutation, f);
    }
  } else {
    f(avail_items, permutation);
    permutation++;
  }
}

template <typename A>
void permute(A items,
             std::function<void(const A &, int)> &f) {
  int permutation = 0;
  permute_helper(items, A::size() - 1, permutation, f);
}

TEST_CASE("Unit Right Triangle", "[HOT]") {
  constexpr const int num_permutations =
      CTMath::partialFactorial(1, tri_verts);
  constexpr const int num_shifts = 10;
  constexpr const int num_tris =
      num_shifts * num_permutations;
  // TODO: Reduce the large relative error of the
  // Wasserstein calculation
  // 19 bits of double precision's 53 is too much!
  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 524288.0;

  std::random_device rd;
  RNG engine(rd());
  constexpr const double min_shift = -10.0;
  constexpr const double max_shift = 9.0;
  // Use float to reduce the number of bits required to
  // represent the final number
  std::uniform_real_distribution<float> genShift(min_shift,
                                                 max_shift);

  Array<Point, 3> points = {Point(0, 0), Point(1, 0),
                            Point(0, 1)};
  Array<Triangle, num_tris> tri;
  std::function<void(const Array<Point, 3> &, int)> filler(
      [&](const Array<Point, 3> &pt_permute,
          int permutation) {
        for(int i = 0; i < num_shifts; i++) {
          double dx = genShift(engine);
          double dy = genShift(engine);
          Array<Point, 3> test_point = {
              Point(pt_permute[0][0] + dx,
                    pt_permute[0][1] + dy),
              Point(pt_permute[1][0] + dx,
                    pt_permute[1][1] + dy),
              Point(pt_permute[2][0] + dx,
                    pt_permute[2][1] + dy)};
          tri[i * num_permutations + permutation] =
              Triangle(test_point[0], test_point[1],
                       test_point[2]);
        }
      });
  permute(points, filler);

  SECTION("Split") {
    for(int i = 0; i < num_tris; i++) {
      boost::variant<std::array<Triangle, 2>,
                     std::array<Triangle, 1> >
          bounds = integral_bounds(tri[i]);
      REQUIRE(bounds.which() == 1);
    }
  }
  SECTION("Wasserstein 2") {
    for(int i = 0; i < num_tris; i++) {
      K_real energy = triangle_w<2>(tri[i]);

      REQUIRE(std::abs(energy * 18 - 1) < max_rel_error);
    }
  }
}

TEST_CASE("Unit Equilateral Triangle", "[HOT]") {
  constexpr const int num_permutations =
      CTMath::partialFactorial(1, tri_verts);
  constexpr const int num_shifts = 10;
  constexpr const int num_tris =
      num_shifts * num_permutations;
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
  std::uniform_real_distribution<float> genShift(min_shift,
                                                 max_shift);

  Array<Point, 3> points = {Point(0, 0), Point(1, 0),
                            Point(0.5, std::sqrt(3) / 2.0)};
  Array<Triangle, num_tris> tri;
  std::function<void(const Array<Point, 3> &, int)> filler(
      [&](const Array<Point, 3> &pt_permute,
          int permutation) {
        for(int i = 0; i < num_shifts; i++) {
          double dx = genShift(engine);
          double dy = genShift(engine);
          Array<Point, 3> test_point = {
              Point(pt_permute[0][0] + dx,
                    pt_permute[0][1] + dy),
              Point(pt_permute[1][0] + dx,
                    pt_permute[1][1] + dy),
              Point(pt_permute[2][0] + dx,
                    pt_permute[2][1] + dy)};
          tri[i * num_permutations + permutation] =
              Triangle(test_point[0], test_point[1],
                       test_point[2]);
        }
      });
  permute(points, filler);

  SECTION("Split") {
    for(int i = 0; i < num_tris; i++) {
      boost::variant<std::array<Triangle, 2>,
                     std::array<Triangle, 1> >
          bounds = integral_bounds(tri[i]);
      REQUIRE(bounds.which() == 0);
    }
  }
  SECTION("Wasserstein 2") {
    for(int i = 0; i < num_tris; i++) {
      K_real energy = triangle_w<2>(tri[i]);
      REQUIRE(std::abs(energy - 0.036084391824351578232) <
              0.036084391824351578232 * max_rel_error);
    }
  }
}
