
#include <iostream>

#include <hot.hpp>

#include "catch.hpp"

TEST_CASE("Triangle Split", "[HOT]") {
  Triangle tri(Point(2, 1), Point(3, 1), Point(2, 2));
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  REQUIRE(bounds.which() == 1);
  auto area = boost::get<std::array<Triangle, 1> >(bounds);
  K_real energy = triangle_w<2>(tri);
  REQUIRE(energy == 0.5 / 9.0);
}
