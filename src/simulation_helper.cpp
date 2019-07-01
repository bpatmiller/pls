#include "simulation.h"

void Simulation::add_sphere_phi() {
  // create a sphere with radius of 1/4 container size
  glm::vec3 center(0.5f * nx * h, 0.5f * ny * h, 0.5f * nz * h);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
        liquid_phi(i, j, k) =
            (glm::distance(center, position) / lx) * 2.0f - 0.5f;
      }
    }
  }
}

void Simulation::intialize_boundaries() {
  for (int i = 0; i < solid_phi.sx; i++) {
    for (int k = 0; k < solid_phi.sy; k++) {
      solid_phi(i, solid_phi.sy - 1, k) = -0.5f * h;
      solid_phi(i, 0, k) = -0.5f * h;
    }
  }

  for (int j = 0; j < solid_phi.sy; j++) {
    for (int k = 0; k < solid_phi.sz; k++) {
      solid_phi(solid_phi.sx - 1, j, k) = -0.5f * h;
      solid_phi(0, j, k) = -0.5f * h;
    }
  }

  for (int i = 0; i < solid_phi.sx; i++) {
    for (int j = 0; j < solid_phi.sy; j++) {
      solid_phi(i, j, solid_phi.sz - 1) = -0.5f * h;
      solid_phi(i, j, 0) = -0.5f * h;
    }
  }
}

// given a position, return the trilinear interpolation
// of the velocity field at that position
glm::vec3 Simulation::trilerp_uvw(glm::vec3 p) {
  glm::ivec3 index;
  glm::vec3 coords;
  glm::vec3 result;
  // u
  position_to_grid(p, U_OFFSET, index, coords);
  result.x = u.trilerp(index, coords);
  // v
  position_to_grid(p, V_OFFSET, index, coords);
  result.y = v.trilerp(index, coords);
  // w
  position_to_grid(p, W_OFFSET, index, coords);
  result.z = w.trilerp(index, coords);
  return result;
}

// does the same as above, but rounds down for easier grid transfer
// use offset = glm::vec3(0) for values sampled at center
void Simulation::position_to_grid(glm::vec3 p, glm::vec3 offset,
                                  glm::ivec3 &index, glm::vec3 &coords) {
  float nx = (p.x / h) - 0.5f + offset.x;
  float ny = (p.y / h) - 0.5f + offset.y;
  float nz = (p.z / h) - 0.5f + offset.z;

  int i = static_cast<int>(nx);
  int j = static_cast<int>(ny);
  int k = static_cast<int>(nz);

  // set index
  index = glm::ivec3(i, j, k);

  float bx = nx - std::floor(nx);
  float by = ny - std::floor(ny);
  float bz = nz - std::floor(nz);

  // set barycentric coordinates
  coords = glm::vec3(bx, by, bz);
}
