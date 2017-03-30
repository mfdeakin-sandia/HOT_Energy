
#include <hot.hpp>

void wasserstein_edge_edge_test() {
  // create two test edges
  Point p0(0, 0);
  Point p1(2, 0);
  Point p2(1, 0);
  Point p3(1, 2);

  // use P10--P11 for both, so distance should be zero
  Segment e0(p0, p1);
  Segment e1(p0, p1);

  // compute distance
  //  wasserstein2_edge_edge( Segment(p0,p1), Segment(p0,p1)
  //  );
  //  wasserstein2_edge_edge( Segment(p0,p1), Segment(p1,p0)
  //  );

  // perp edges at midpoint
  // of uniform length
  // appears that W2 is invariant to permutation of matching
  // up integration points!
  // W1 order matters, best is a large spread, match close
  // points together and match far points together
  // W1 best is interlaced, alternating points from each
  // side, *not* 0-5 matched with 0-5
  if(0) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, 1), Point(1, 3)));  // 2.22  W2
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, 0), Point(1, 2)));  // 1.39
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -0.3), Point(1, 1.7)));  // 1.19
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -0.5), Point(1, 1.5)));  // 1.09
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(1, -1), Point(1, 1)));  // 0.966
  }

  // perp edges, not at midpoint
  // appears that W2 is invariant to permutation of matching
  // up integration points!
  // W1 is as in the above case
  if(0) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, 1), Point(0, 3)));  // 2.43 W2
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, 0), Point(0, 2)));  // 1.71
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -0.3), Point(0, 1.7)));  // 1.55
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -0.5), Point(0, 1.5)));  // 1.47
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0, -1), Point(0, 1)));  // 1.39
  }

  // skew
  if(1) {
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 1)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 2)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(-1, 1), Point(1, 3)));

    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 1)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 2)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, 1), Point(1.8, 3)));

    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 0.3)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 0.9)));
    wasserstein2_edge_edge(
        Segment(p0, p1),
        Segment(Point(0.2, -0.3), Point(1.8, 1.5)));
  }
}

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
  // Test Wasserstein distance between two segments
  if(0) {
    wasserstein_edge_edge_test();
    test_tri_w2();
    return 0;
  }

  // Test hot energy between triangle and triangle* == point

  DT dt;
  const int num_points = 30;
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
  return 0;
}
