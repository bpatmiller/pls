#pragma once

#include "aabb_triangle_overlap.h"
#include "grid.h"
#include <fstream>
#include <iostream>

// filename of .ply file
// grid to export voxels to
// offset of the file in grid space
void voxelize_mesh(std::string fname, Grid &grid, glm::vec3 offset,
                   int cell_type) {
  std::vector<glm::vec3> vertices;
  std::vector<glm::uvec3> indices;

  std::ifstream in;
  std::string ln;
  int vertex_num = 0;
  int tri_num = 0;

  std::cout << "   loading mesh...";
  in.open(fname);

  // parse header
  while (std::getline(in, ln)) {
    std::stringstream ss(ln);
    std::string first;
    ss >> first;
    if (first == "element") {
      int n;
      ss >> first >> n;
      if (!vertex_num) {
        vertex_num = n;
      } else {
        tri_num = n;
      }
    } else if (first == "end_header") {
      break;
    }
  }

  // parse vertices
  for (int i = 0; i < vertex_num; i++) {
    float x, y, z;
    std::getline(in, ln);
    std::stringstream ss(ln);
    ss >> x >> y >> z;
    vertices.push_back(glm::vec3(x, y, z) + offset);
  }

  // parse indices
  for (int i = 0; i < tri_num; i++) {
    uint t, a, b, c;
    std::getline(in, ln);
    std::stringstream ss(ln);
    ss >> t >> a >> b >> c;
    indices.push_back(glm::uvec3(a, b, c));
  }

  std::cout << "done (" << vertex_num << " vertices and " << tri_num
            << " triangles)\n   voxelizing mesh...";

  assert((sizeof(glm::vec3)) == 3 * sizeof(float));

  // for each grid node, do AABB intersection
  for (int i = 0; i < grid.marker.sx; i++) {
    for (int j = 0; j < grid.marker.sy; j++) {
      for (int k = 0; k < grid.marker.sz; k++) {
        // brute force check all triangles i guess
        for (int t = 0, tmax = indices.size(); t < tmax; t++) {
          glm::uvec3 ind = indices[t];
          glm::vec3 v1 = vertices[ind[0]];
          glm::vec3 v2 = vertices[ind[1]];
          glm::vec3 v3 = vertices[ind[2]];

          float boxcenter[3] = {(i + 0.5f) * grid.h, (j + 0.5f) * grid.h,
                                (k + 0.5f) * grid.h};
          float boxhalfsize[3] = {grid.h * 0.5f, grid.h * 0.5f, grid.h * 0.5f};
          float triverts[3][3] = {
              {v1.x, v1.y, v1.z}, {v2.x, v2.y, v2.z}, {v3.x, v3.y, v3.z}};

          if (triBoxOverlap(boxcenter, boxhalfsize, triverts)) {
            grid.marker(i, j, k) = cell_type;
          }
        }
      }
    }
  }

  std::cout << "done\n";
}