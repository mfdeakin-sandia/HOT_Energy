// ply_writer.cpp
#include "ply_writer.hpp"

#include <iostream>
#include <fstream>

void aggregate_simplices(const DT &mesh, std::list<Point> &points,
                         std::list<std::array<int, tri_verts> > &faces) {
  std::map<Point, int> point_indices;
  for (auto face_itr = mesh.finite_faces_begin();
       face_itr != mesh.finite_faces_end(); face_itr++) {
    std::array<int, tri_verts> face_verts;

    for (int i = 0; i < tri_verts; i++) {
      Point p = face_itr->vertex(i)->point();
      auto key_itr = point_indices.find(p);
      if (key_itr != point_indices.end()) {
        face_verts[i] = key_itr->second;
      } else {
        int point_idx = points.size();
        points.push_back(p);
        point_indices[p] = point_idx;
        face_verts[i] = point_idx;
      }
    }
    faces.push_front(face_verts);
  }
}

void write_ply(const char *fname, const DT &mesh) {
  std::ofstream output(fname);
  static constexpr const char *header = "ply\n"
                                        "format ascii 1.0";
  output << header << std::endl;

  //
  std::list<Point> points;
  std::list<std::array<int, tri_verts> > faces;
  aggregate_simplices(mesh, points, faces);

  static constexpr const char *vertex_props = "property float x\n"
                                              "property float y\n"
                                              "property float z";
  output << "element vertex " << points.size() << std::endl << vertex_props
         << std::endl;

  static constexpr const char *face_props =
      "property list uchar int vertex_index";
  output << "element face " << faces.size() << std::endl << face_props
         << std::endl;

  output << "end_header" << std::endl;

  // Output the vertices
  for(const Point &p : points) {
    output << p << ' ' << 0 << std::endl;
  }
  for (const std::array<int, tri_verts> &f : faces) {
    output << tri_verts;
    for (int point_idx : f) {
      output << " " << point_idx;
    }
    output << std::endl;
  }
}
