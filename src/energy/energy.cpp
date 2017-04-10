
#include <hot.hpp>
#include <ply_writer.hpp>

#include <random>

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

void generate_rand_dt(int num_points, DT &dt) {
  std::random_device rd;
  RNG engine(rd());
  std::uniform_real_distribution<double> genPos(min_pos,
                                                max_pos);
  for(int i = 0; i < num_points; i++) {
    Point pt(genPos(engine), genPos(engine));
    dt.insert(pt);
  }
}

int main(int argc, char **argv) {

  // Test hot energy between triangle and triangle* == point

  DT dt;
  const int num_points = 10;
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
  K_real energy = hot_energy<2>(dt);
  std::cout << "Mesh energy: " << energy << std::endl;
  write_ply("test.ply", dt);
  hot_optimize<2>(dt);
  return 0;
}
