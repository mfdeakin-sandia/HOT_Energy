
#include <cmath>
#include <limits>

#include <iomanip>
#include <iostream>

#include <hot.hpp>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Right Triangle", "[HOT]") {
  // TODO: Reduce the large relative error of this
  // calculation
  // 10 bits of double precision is too much!
  constexpr const double max_rel_error =
      std::numeric_limits<double>::epsilon() * 1024.0;

  Triangle tri(Point(2, 1), Point(3, 1), Point(2, 2));
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  SECTION("Split") { REQUIRE(bounds.which() == 1); }
  SECTION("HOT Energy") {
    auto area =
        boost::get<std::array<Triangle, 1> >(bounds);
    K_real energy = triangle_w<2>(tri);

    REQUIRE(std::abs(energy * 18 - 1) < 18 * max_rel_error);
  }
}
