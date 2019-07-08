#include "simulation.h"

void Simulation::add_fluid(float density_) {
  fluids.push_back(Fluid(density_, nx, ny, nz));
}

void Simulation::add_pool_phi() {
  fluids.clear();
  add_fluid(1.0f);
  add_fluid(0.0f);
  // 0 density implies we do not include this in our computations
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        float phi_val = (h * j) - (0.5f * ly);
        fluids[0].phi(i, j, k) = phi_val;
        fluids[1].phi(i, j, k) = -phi_val;
      }
    }
  }
  enforce_boundaries();
}

void Simulation::add_pool_drop() {
  fluids.clear();
  add_fluid(1.0f);
  add_fluid(2.0f);
  add_fluid(0.0f);

  glm::vec3 center(0.5f * (nx)*h, 0.7f * (ny)*h, 0.5f * (nz)*h);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
        float phi_val_pool = (h * j) - (0.25f * ly);
        float phi_val_sphere = glm::distance(center, position) - (lx * 0.2f);
        fluids[0].phi(i, j, k) = phi_val_sphere;
        fluids[1].phi(i, j, k) = phi_val_pool;
        fluids[2].phi(i, j, k) = std::max(-phi_val_pool, -phi_val_sphere);
      }
    }
  }
  enforce_boundaries();
}

void Simulation::add_two_walls() {
  fluids.clear();
  add_fluid(1.0f);
  add_fluid(2.0f);
  add_fluid(0.0f);
  // 0 density implies we do not include this in our computations
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        float phi_val_l = (h * i) - (0.4f * lx);
        float phi_val_r = (lx - h * i) - (0.4f * lx);

        fluids[0].phi(i, j, k) = phi_val_l;
        fluids[1].phi(i, j, k) = phi_val_r;
        fluids[2].phi(i, j, k) = std::max(-phi_val_l, -phi_val_r);
      }
    }
  }
  enforce_boundaries();
}

void Simulation::add_two_spheres() {
  fluids.clear();
  add_fluid(1.0f);
  add_fluid(2.0f);
  add_fluid(0.0f);
  glm::vec3 center_l(0.3f * (nx)*h, 0.7f * (ny)*h, 0.5f * (nz)*h);
  glm::vec3 center_r(0.7f * (nx)*h, 0.7f * (ny)*h, 0.5f * (nz)*h);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
        float phi_val_l = glm::distance(center_l, position) - (lx * 0.2f);
        float phi_val_r = glm::distance(center_r, position) - (lx * 0.2f);

        fluids[0].phi(i, j, k) = phi_val_l;
        fluids[1].phi(i, j, k) = phi_val_r;
        fluids[2].phi(i, j, k) = std::max(-phi_val_l, -phi_val_r);
      }
    }
  }
}

void Simulation::add_sphere_phi() {
  fluids.clear();
  add_fluid(1.0f);
  add_fluid(0.0f);
  // create a sphere with radius of 1/4 container size
  glm::vec3 center(0.5f * (nx)*h, 0.5f * (ny)*h, 0.5f * (nz)*h);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
        float phi_val = glm::distance(center, position) - (lx * 0.3f);
        fluids[0].phi(i, j, k) = phi_val;
        fluids[1].phi(i, j, k) = -phi_val;
      }
    }
  }
}

// inigo quile'z capsule distance function
static inline float sdCapsule(glm::vec3 p, glm::vec3 a, glm::vec3 b, float r) {
  glm::vec3 pa = p - a, ba = b - a;
  float h = glm::clamp(glm::dot(pa, ba) / glm::dot(ba, ba), 0.0f, 1.0f);
  return glm::length(pa - ba * h) - r;
}

// void Simulation::add_cylinder_phi() {
//   glm::vec3 a(h * (nx * 0.35f + 0.5f), h * (nx * 0.5f + 0.5f),
//               h * (nz * 0.5f + 0.5f));
//   glm::vec3 b(h * (nx * 0.65f + 0.5f), h * (nx * 0.5f + 0.5f),
//               h * (nz * 0.5f + 0.5f));
//   float r = lx * 0.25f;
//   for (int i = 0; i < nx; i++) {
//     for (int j = 0; j < ny; j++) {
//       for (int k = 0; k < nz; k++) {
//         glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);
//         liquid_phi(i, j, k) = sdCapsule(position, a, b, r);
//       }
//     }
//   }
// }

void Simulation::set_vortex_field() {
  // predefined_field = true;
  v.set(0);
  glm::vec2 center(lx * 0.5f, ly * 0.5f);
  for (int i = 0; i < nx + 1; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz + 1; k++) {
        glm::vec2 position((i + 0.5f) * h, (k + 0.5f) * h);
        float x = position.x - center.x;
        float y = position.y - center.y;
        if (k < nz)
          u(i, j, k) = 5.0 * -y / (x * x + y * y);
        if (i < nx)
          w(i, j, k) = 5.0 * x / (x * x + y * y);
      }
    }
  }
}

void Simulation::intialize_boundaries() {
  solid_phi.set(0.5f * h);
  for (int i = 0; i < solid_phi.sx; i++) {
    for (int k = 0; k < solid_phi.sy; k++) {
      solid_phi(i, solid_phi.sy - 1, k) = -0.5f * h;
      solid_phi(i, solid_phi.sy - 2, k) = -0.5f * h;
      solid_phi(i, 0, k) = -0.5f * h;
      solid_phi(i, 1, k) = -0.5f * h;
    }
  }

  for (int j = 0; j < solid_phi.sy; j++) {
    for (int k = 0; k < solid_phi.sz; k++) {
      solid_phi(solid_phi.sx - 1, j, k) = -0.5f * h;
      solid_phi(solid_phi.sx - 2, j, k) = -0.5f * h;
      solid_phi(0, j, k) = -0.5f * h;
      solid_phi(1, j, k) = -0.5f * h;
    }
  }

  for (int i = 0; i < solid_phi.sx; i++) {
    for (int j = 0; j < solid_phi.sy; j++) {
      solid_phi(i, j, solid_phi.sz - 1) = -0.5f * h;
      solid_phi(i, j, solid_phi.sz - 2) = -0.5f * h;
      solid_phi(i, j, 0) = -0.5f * h;
      solid_phi(i, j, 1) = -0.5f * h;
    }
  }

  Fluid f(0.0, nx, ny, nz);
  f.phi = solid_phi;
  reinitialize_phi(f);
  solid_phi = f.phi;
}

// given a position, return the trilinear interpolation
// of the velocity field at that position
glm::vec3 Simulation::trilerp_tutvtw(glm::vec3 p) {
  p = glm::clamp(p, 0.0001f * h, lx - (1.0001f * h));
  glm::ivec3 index;
  glm::vec3 coords;
  glm::vec3 result;
  // u
  position_to_grid(p, U_OFFSET, index, coords);
  result.x = tu.trilerp(index, coords);
  // v
  position_to_grid(p, V_OFFSET, index, coords);
  result.y = tv.trilerp(index, coords);
  // w
  position_to_grid(p, W_OFFSET, index, coords);
  result.z = tw.trilerp(index, coords);
  return result;
}

glm::vec3 Simulation::trilerp_uvw(glm::vec3 p) {
  p = glm::clamp(p, 0.5001f * h, lx - (0.5001f * h));
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

  // p = glm::clamp(p, 0.5001f * h, lx - (0.5001f * h));

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
