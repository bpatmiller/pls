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

// check all grid cells with |phi|<3h, seed 64 particles in each
void Simulation::reseed_particles() {
  // set the valid flag if particles are close to the interface
  for (auto &p : particles) {
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, CENTER_OFFSET, index, coords);
    float phi_val = liquid_phi.trilerp(index, coords);
    if (std::abs(phi_val) < 3.0f * h) {
      p.valid = true;
    } else {
      p.valid = false;
      particle_count(index.x, index.y, index.z) -= 1;
    }
  }
  // delete particles that are far away from the interface
  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [](Particle const &p) { return !p.valid; }),
                  particles.end());

  // make sure each grid cell near the interface has 64 particles
  for (int i = 0; i < liquid_phi.sx; i++) {
    for (int j = 0; j < liquid_phi.sy; j++) {
      for (int k = 0; k < liquid_phi.sz; k++) {
        if (std::abs(liquid_phi(i, j, k)) < 3.0f * h) {
          // reseed
          for (int n = particle_count(i, j, k); n < 64; n++) {
            float jitterx = glm::linearRand(0.0f, 1.0f);
            float jittery = glm::linearRand(0.0f, 1.0f);
            float jitterz = glm::linearRand(0.0f, 1.0f);
            glm::vec3 position =
                glm::vec3((float)i + jitterx, (float)j + jittery,
                          (float)k + jitterz) *
                h;
            Particle p(position, liquid_phi.trilerp(
                                     glm::ivec3(i, j, k),
                                     glm::vec3(jitterx, jittery, jitterz)));
            particles.push_back(p);
            particle_count(i, j, k) += 1;
          }
        }
      }
    }
  }
}

void Simulation::compute_derivatives() {
  // compute foward, backward, and central difference derivatives of phi
  for (int i = 0; i < liquid_phi.sx; i++) {
    for (int j = 0; j < liquid_phi.sy; j++) {
      for (int k = 0; k < liquid_phi.sz; k++) {
        // dx forward, backward, center
        if (i < liquid_phi.sx - 1)
          DX(i, j, k).x = (liquid_phi(i + 1, j, k) - liquid_phi(i, j, k)) / h;
        if (i > 0)
          DX(i, j, k).y = (liquid_phi(i, j, k) - liquid_phi(i - 1, j, k)) / h;
        if (i < liquid_phi.sx - 1 and i > 0)
          DX(i, j, k).z =
              (liquid_phi(i + 1, j, k) - liquid_phi(i - 1, j, k)) / (2.0f * h);
        // dy forward, backward, center
        if (j < liquid_phi.sy - 1)
          DY(i, j, k).x = (liquid_phi(i, j + 1, k) - liquid_phi(i, j, k)) / h;
        if (j > 0)
          DY(i, j, k).y = (liquid_phi(i, j, k) - liquid_phi(i, j - 1, k)) / h;
        if (j < liquid_phi.sy - 1 and j > 0)
          DY(i, j, k).z =
              (liquid_phi(i, j + 1, k) - liquid_phi(i, j - 1, k)) / (2.0f * h);
        // dz forward, backward, center
        if (k < liquid_phi.sz - 1)
          DZ(i, j, k).x = (liquid_phi(i, j, k + 1) - liquid_phi(i, j, k)) / h;
        if (k > 0)
          DZ(i, j, k).y = (liquid_phi(i, j, k) - liquid_phi(i, j, k - 1)) / h;
        if (k < liquid_phi.sz - 1 and k > 0)
          DZ(i, j, k).z =
              (liquid_phi(i, j, k + 1) - liquid_phi(i, j, k - 1)) / (2.0f * h);
      }
    }
  }
}

// compute gradient norm via Godunov scheme
void Simulation::norm_gradient() {
  norm_grad.clear();
  for (int i = 0; i < norm_grad.sx; i++) {
    for (int j = 0; j < norm_grad.sy; j++) {
      for (int k = 0; k < norm_grad.sz; k++) {
      }
    }
  }
}

void Simulation::reinitialize_phi() {
  // compute sigmoid function for phi_0
  sig = liquid_phi / (liquid_phi * liquid_phi + h * h);

  norm_gradient();
}

// particle advection with bounds checking
void Simulation::advect_particles(float dt) {
  for (auto &p : particles) {
    p.position = rk2(p.position, dt);
    // bounds checking
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, CENTER_OFFSET, index, coords);
    float solid_phi_val = solid_phi.trilerp(index, coords);
    if (false && solid_phi_val < 0) {
      glm::vec3 grad = solid_phi.gradlerp(index, coords, h);
      if (glm::length(grad) > 0) {
        grad = glm::normalize(grad);
        p.position -= solid_phi_val * grad;
      }
    }
  }
}

// second-order accurate runge-kutta advection
glm::vec3 Simulation::rk2(glm::vec3 position, float dt) {
  float min_x = h * 1.0001f;
  float max_x = h * ((float)nx - 1.0001f);

  // first stage of RK2
  glm::vec3 gu = trilerp_uvw(position);
  glm::vec3 mid = position + 0.5f * dt * gu;
  mid = glm::clamp(mid, min_x, max_x);
  // second stage
  gu = trilerp_uvw(mid);
  mid = position + dt * gu;
  mid = glm::clamp(mid, min_x, max_x);
  return mid;
}

// advect phi with equation phi_t = - V dot grad phi
void Simulation::advect_phi(float dt) {
  phi_copy.clear();
  for (int i = 1; i < liquid_phi.sx - 1; i++) {
    for (int j = 1; j < liquid_phi.sy - 1; j++) {
      for (int k = 1; k < liquid_phi.sz - 1; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);

        glm::vec3 velocity = trilerp_uvw(position);
        glm::vec3 phi_grad = liquid_phi.upwind_gradient(i, j, k, h, velocity);

        phi_copy(i, j, k) =
            liquid_phi(i, j, k) - dt * glm::dot(velocity, phi_grad);
      }
    }
  }
  phi_copy.copy_to(liquid_phi);
}

// semi-lagrangian advection via backwards euler as in [Stam 1999]
void Simulation::advect_velocity(float dt) {
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
void Simulation::add_gravity(float dt) { v -= 9.8 * dt; }

void Simulation::enforce_boundaries() {
  for (int i = 0; i < solid_phi.sx; i++) {
    for (int j = 0; j < solid_phi.sy; j++) {
      for (int k = 0; k < solid_phi.sz; k++) {
        // if solid cell, set face velocities to 0
        if (solid_phi(i, j, k) < 0) {
          liquid_phi(i, j, k) = 0.5f * h;
          u(i, j, k) = 0;
          u(i + 1, j, k) = 0;
          v(i, j, k) = 0;
          v(i, j + 1, k) = 0;
          w(i, j, k) = 0;
          w(i, j, k + 1) = 0;
        }
      }
    }
  }
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
    // reseed_particles();
    reinitialize_phi();
    advect_particles(substep);
    advect_phi(substep);
    advect_velocity(substep);
    add_gravity(substep);
    enforce_boundaries();
    project(substep);
    // extrapolate(u, u_valid);
    // extrapolate(v, v_valid);
    // extrapolate(w, w_valid);
    // enforce_boundaries();
  }
}