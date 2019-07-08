#pragma once

#include "aabb_triangle_overlap.h"
#include <fstream>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>

class Mesh {
public:
  std::vector<glm::vec3> vertices;
  std::vector<glm::uvec3> indices;
  std::string base_filename;

  Mesh() { base_filename = "filename_"; }

  void save(std::string fname) {
    std::ofstream ofile(std::string("out/") + base_filename + fname);

    ofile << "ply\n";
    ofile << "format ascii 1.0\n";
    ofile << "element vertex " << vertices.size() << "\n";
    ofile << "property float x\n";
    ofile << "property float y\n";
    ofile << "property float z\n";
    ofile << "element face " << indices.size() << "\n";
    ofile << "property list uchar int vertex_index\n";
    ofile << "end_header\n";

    for (auto &v : vertices) {
      ofile << v.x << " " << v.z << " " << v.y << "\n";
    }
    for (auto &i : indices) {
      ofile << "3 " << i.x << " " << i.y << " " << i.z << "\n";
    }

    ofile.close();
    std::cout << "saved out/m_" << fname << " containing " << vertices.size()
              << " vertices\n";
  }
};