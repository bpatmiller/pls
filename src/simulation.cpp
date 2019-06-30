#include "simulation.h"

static inline double sqr(double d) { return std::pow(d, 2.0); }
static inline float sqr(float f) { return std::pow(f, 2.0f); }
static inline float three_max(float a, float b, float c) {
  return std::max(a, std::max(b, c));
}

float Simulation::CFL() {
  float max_v = three_max(u.infnorm(), w.infnorm(), v.infnorm());
  max_v = std::max(max_v, 1e-16f);
  return h / max_v;
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

void Simulation::advect_particles(float dt) {
  glm::vec3 mid, gu;
  // boundaries
  float xmin = 1.01 * h;
  float ymin = 1.01 * h;
  float zmin = 1.01 * h;
  float xmax = lx - 1.01 * h;
  float ymax = ly - 1.01 * h;
  float zmax = lz - 1.01 * h;
  for (uint i = 0; i < particles.size(); i++) {
    Particle &p = particles[i];
    // first stage of RK2
    gu = trilerp_uvw(p.position);
    mid = p.position + 0.5f * dt * gu;
    mid.x = glm::clamp(mid.x, xmin, xmax);
    mid.y = glm::clamp(mid.y, ymin, ymax);
    mid.z = glm::clamp(mid.z, zmin, zmax);
    // second stage
    gu = trilerp_uvw(mid);
    p.position += dt * gu;
    p.position.x = glm::clamp(p.position.x, xmin, xmax);
    p.position.y = glm::clamp(p.position.y, ymin, ymax);
    p.position.z = glm::clamp(p.position.z, zmin, zmax);
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
      solid_phi(i, j, solid_phi.sz - 1) = -0.5f *h;
      solid_phi(i, j, 0) = -0.5f * h;
    }
  }
}

void Simulation::advance(float dt) {
  float t = 0.0f;
  float substep;
  while (t < dt) {
    substep = CFL();
    if (t + substep >= dt) {
      substep = dt - t;
    }
    t += substep;
    // perform each simulation step
    advect_particles(substep);
    // advect(substep);
    // add_gravity(substep);
    // project(substep);
    // extrapolate(u, u_valid);
    // extrapolate(v, v_valid);
    // extrapolate(w, w_valid);
    // enforce_boundaries();
  }
}