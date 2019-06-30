#pragma once
#include "array3.h"
#include <glm/glm.hpp>
#include <vector>

// struct GRIDCELL {
//    glm::vec3 p[8];		//position of each corner of the grid in world
//    space float val[8];	//value of the function at this grid corner
// };

// given a grid cell, returns the set of triangles that approximates the region
// where val == 0.
void polygonize(std::vector<glm::vec3> &positions, std::vector<float> &values,
                int &vert_count, std::vector<glm::uvec3> &indices,
                std::vector<glm::vec3> &vertices);