
#ifndef _OPTIMIZED_MESH_HPP_
#define _OPTIMIZED_MESH_HPP_

#include <limits>

#include <array>
#include <list>

static constexpr const int space_dim = 2;

/* OptimizedMesh implements Newton's method to optimize the mesh according to
 * the specified gradient function
 * OptimizedMesh inherits from it's template paramter Mesh,
 * which must inherit from CGAL::Triangulation
 */
template <typename Mesh, typename real> class OptimizedMesh : public Mesh {
public:
  OptimizedMesh(Triangulation &mesh, real stop_threshold) : Mesh(mesh) {
    optimize(stop_threshold);
  }

  std::list<Mesh::Vertex_handle> internal_vertices() {
    std::list<Mesh::Vertex_handle> verts;
    std::set<Mesh::Vertex_handle> seen;
    // First mark the external vertices as having been seen
    for (auto vert_itr = incident_vertices(infinite_vertex());
         seen.count(vert_itr) == 0; vert_itr++) {
      seen.insert(vert_itr);
    }
    // Now iterator over all of the vertices,
    // storing them if we haven't seen them before
    for (auto vert_itr = finite_vertices_begin();
         vert_itr != dt.finite_vertices_end(); vert_itr++) {
      if (seen.count(vert_itr) == 0) {
        seen.insert(vert_itr);
        verts.push_back(vert_itr);
      }
    }
    return verts;
  }

private:
  virtual void optimize(real stop_threshold) {
    double prev_energy = -std::numeric_limits<real>::infinity();
    double cur_energy = energy();
    const real scale = 1.0;
    while ((cur_energy - prev_energy) > stop_threshold) {
      for (Mesh::Vertex_handle vtx : internal_vertices()) {
        std::array<real, space_dim> gradient = vertex_energy_grad(vtx);
        Mesh::Point updated_point = vtx->point();
        for (int i = 0; i < space_dim; i++) {
          updated_point[i] -= scale * gradient[i];
        }
        updated_point move(vtx, updated_point);
      }
      prev_energy = cur_energy;
      cur_energy = energy();
    }
  }

  virtual real energy() = 0;
  virtual std::array<real, space_dim>
  vertex_energy_grad(Mesh::Vertex_handle vtx) = 0;

  static Triangle face_to_tri(const Face &face) {
    return Triangle(face.vertex(0)->point(), face.vertex(1)->point(),
                    face.vertex(2)->point());
  }
};

#endif
