#include <hot.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


#define PI 3.14159265

void test_tri_w2() {
  Triangle tri(Point(2, 1), Point(3, 1), Point(2, 2));
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  assert(bounds.which() == 1);
  auto area = boost::get<std::array<Triangle, 1> >(bounds);

  std::cout
      << std::endl
      << "Computing Wasserstein distance to the dual point"
      << std::endl
      << std::endl;
  K_real value = triangle_w<2>(tri);
  std::cout << "Unit Right Triangle expected Wasserstein "
               "Distance: 5/90 (~0.0555555). Calculated "
               "Distance: "
            << value << std::endl
            << std::endl;
}

using RNG = std::mt19937_64;

constexpr const float min_pos = 10.0;
constexpr const float max_pos = 20.0;

int main(int argc, char **argv) {
	Triangle tri(Point(0,0), Point(1,0), Point(0,1));
	double star0_energy=tri_energy<2,0>(tri); 
	double star1_energy=tri_energy<2,1>(tri);
	double star2_energy=tri_energy<2,2>(tri);

	K_real tri_was=triangle_w<2>(tri);
	double area=tri.area();
	Mich_energy=area*tri_was; 
		
	std::cout<< std::setw(15) << star_2_energy <<std::setw(15) << Mich_energy <<std::endl; 

	return 0;
}
