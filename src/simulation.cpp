#include "simulation.h"

// ----------------
// simulation steps
// ----------------

// returns a safe timestep
float Simulation::CFL() {
  float max_v = std::max(u.infnorm(), std::max(w.infnorm(), v.infnorm()));
  max_v = std::max(max_v, 1e-16f);
  return h / max_v;
}

// particle advection with bounds checking
void Simulation::advect_particles(float dt) {
  for (uint i = 0; i < particles.size(); i++) {
    Particle &p = particles[i];
    p.position = rk2(p.position, dt);
    // bounds checking
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, CENTER_OFFSET, index, coords);
    float solid_phi_val = solid_phi.trilerp(index, coords);
    if (solid_phi_val < 0) {
      glm::vec3 grad = solid_phi.gradlerp(index, coords);
      if (glm::length(grad) > 0) {
        grad = glm::normalize(grad);
        p.position -= solid_phi_val * grad;
      }
    }
  }
}

// second-order accurate runge-kutta advection
glm::vec3 Simulation::rk2(glm::vec3 position, float dt) {
  // first stage of RK2
  glm::vec3 gu = trilerp_uvw(position);
  glm::vec3 mid = position + 0.5f * dt * gu;
  // second stage
  gu = trilerp_uvw(mid);
  return position + dt * gu;
}

// semi-lagrangian advection via backwards euler
void Simulation::advect(float dt) {
  tu.clear();
  tv.clear();
  tw.clear();

  // advect u
  for (int i = 0; i < u.sx; i++) {
    for (int j = 0; j < u.sy; j++) {
      for (int k = 0; k < u.sz; k++) {
        glm::vec3 position(i * h, (j + 0.5f) * h, (k + 0.5f) * h);
        position = rk2(position, -dt);
        tu(i, j, k) = trilerp_uvw(position).x;
      }
    }
  }

  // advect v
  for (int i = 0; i < v.sx; i++) {
    for (int j = 0; j < v.sy; j++) {
      for (int k = 0; k < v.sz; k++) {
        glm::vec3 position((i + 0.5f) * h, j * h, (k + 0.5f) * h);
        position = rk2(position, -dt);
        tv(i, j, k) = trilerp_uvw(position).y;
      }
    }
  }

  // advect w
  for (int i = 0; i < w.sx; i++) {
    for (int j = 0; j < w.sy; j++) {
      for (int k = 0; k < w.sz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, k * h);
        position = rk2(position, -dt);
        tw(i, j, k) = trilerp_uvw(position).z;
      }
    }
  }

  tu.copy_to(u);
  tv.copy_to(v);
  tw.copy_to(w);
}

// apply acceleration due to gravity
void Simulation::add_gravity(float dt) {
  for (int i = 0; i < v.size; i++) {
    v.data[i] -= 9.8f * dt;
  }
}

// compute signed distance function from fluid
void Simulation::compute_phi() {
  liquid_phi.set(3 * h);
  for (auto &p : particles) {
    glm::ivec3 index;
    glm::vec3 coords;
    position_to_grid(p.position, CENTER_OFFSET, index, coords);
    // loop through neighboring cells
    for (int i = std::max((float)0, (float)index.x - 1);
         i <= std::min((float)index.x + 1, (float)nx - 1); i++) {
      for (int j = std::max((float)0, (float)index.y - 1);
           j <= std::min((float)index.y + 1, (float)ny - 1); j++) {
        for (int k = std::max((float)0, (float)index.z - 1);
             k <= std::min((float)index.z + 1, (float)nz - 1); k++) {
          glm::vec3 grid_position((index.x + i + 0.5f) * h,
                                  (index.y + j + 0.5f) * h,
                                  (index.z + k + 0.5f) * h);
          float phi_val =
              glm::distance(p.position, grid_position) - particle_radius;
          if (phi_val < liquid_phi(index.x + i, index.y + j, index.z + k)) {
            liquid_phi(index.x + i, index.y + j, index.z + k) = phi_val;
          }
        }
      }
    }
  }
  // extend phi into solid cells
  // TODO
}

// project the velocity field onto its divergence-free part
void Simulation::project(float dt) {}

// advance the simulation by a given time
// the time will be split up into substeps
// to satisfy the CFL condition
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
    advect(substep);
    add_gravity(substep);
    compute_phi();
    project(substep);
    // extrapolate(u, u_valid);
    // extrapolate(v, v_valid);
    // extrapolate(w, w_valid);
    // enforce_boundaries();
  }
}