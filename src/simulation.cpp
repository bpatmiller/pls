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

float Simulation::trilerp_scalar_field(Array3f &field, glm::vec3 position) {
  glm::vec3 coords;
  glm::ivec3 index;
  position_to_grid(position, glm::vec3(0), index, coords);
  return field.trilerp(index, coords);
}

// check all eight corners of a grid cell for a small phi value
bool Simulation::check_corners(int i, int j, int k) {
  for (int x = 0; x <= 1; x++) {
    for (int y = 0; y <= 1; y++) {
      for (int z = 0; z <= 1; z++) {
        float phi_val = liquid_phi(i + x, j + y, k + z);
        if (std::abs(phi_val) < 3.0f * h)
          return true;
      }
    }
  }
  return false;
}

void Simulation::reseed_particles() {
  if (reseed_counter % 10 == 0) {
    reseed_counter = 1;
    remove_particles();
    initialize_particles();
    // std::cout << "reseeding particles\n";
  }
  reseed_counter++;
}

void Simulation::initialize_particles() {
  int particles_added = 0;
  // count how many particles are in each cell
  particle_count.clear();
  for (auto &p : particles) {
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, glm::vec3(0), index, coords);
    particle_count(index) += 1;
  }

  // "top off" remaining cells with small phi values
  for (int i = 1; i < liquid_phi.sx - 1; i++) {
    for (int j = 1; j < liquid_phi.sy - 1; j++) {
      for (int k = 1; k < liquid_phi.sz - 1; k++) {
        if (check_corners(i, j, k)) {
          glm::vec3 base_position((i + 0.5f) * h, (j + 0.5f) * h,
                                  (k + 0.5f) * h);
          for (int count = particle_count(i, j, k); count < 64; count++) {
            glm::vec3 position =
                base_position + glm::linearRand(glm::vec3(0.0f), glm::vec3(h));
            float phi_val = trilerp_scalar_field(liquid_phi, position);
            // add particle with position, starting phi, and radius
            Particle p(position, phi_val,
                       glm::clamp(std::abs(phi_val), 0.1f * h, 0.5f * h));
            particles.push_back(p);
            particle_count(i, j, k) += 1;
            particles_added++;
          }
        }
      }
    }
  }
  // std::printf("added %i particles\n", particles_added);
}

// remove all extraneous particles
void Simulation::remove_particles() {
  // int old_size = particles.size();
  // mark all distant particles to be deleted
  for (auto &p : particles) {
    float phi_val = trilerp_scalar_field(liquid_phi, p.position);
    p.valid = (std::abs(phi_val) < 3.0f * h);
  }
  // remove from particle vector
  particles.erase(std::remove_if(particles.begin(), particles.end(),
                                 [](Particle const &p) { return !p.valid; }),
                  particles.end());

  // std::printf("removed %i particles\n", old_size - (int)particles.size());
  // std::printf("particle.size()=%i\n", (int)particles.size());
}

void Simulation::correct_levelset() {
  // FIXME check equals operator
  liquid_phi_minus = liquid_phi;
  liquid_phi_plus = liquid_phi;
  for (auto &p : particles) {
    // determine if the particle is escaped
    float local_phi_val = trilerp_scalar_field(liquid_phi, p.position);
    if (p.starting_phi * local_phi_val < 0 && local_phi_val > p.radius) {
      // compute the phi_p at all eight grid points
      glm::ivec3 index;
      glm::vec3 coords;
      position_to_grid(p.position, glm::vec3(0), index, coords);
      float sign_p = (p.starting_phi > 0) ? 1.0f : -1.0f;
      // compare phi_p to each neighboring grid node
      for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 1; j++) {
          for (int k = 0; k <= 1; k++) {
            glm::vec3 grid_position((index.x + i + 0.5f) * h,
                                    (index.y + j + 0.5f) * h,
                                    (index.z + k + 0.5f) * h);
            // equation 10 in PLS
            float phi_p =
                sign_p * (p.radius - glm::distance(grid_position, p.position));
            // rebuild phi_+
            if (sign_p > 0) {
              liquid_phi_plus(index.x + i, index.y + j, index.z + k) =
                  std::max(phi_p, liquid_phi_plus(index.x + i, index.y + j,
                                                  index.z + k));
            }
            // rebuild phi_-
            else {
              liquid_phi_plus(index.x + i, index.y + j, index.z + k) =
                  std::min(phi_p, liquid_phi_plus(index.x + i, index.y + j,
                                                  index.z + k));
            }
          }
        }
      }
    }
  }

  // merge phi_+ and phi_-
  // equation 13 in PLS
  for (int i = 0; i < liquid_phi.size; i++) {
    if (std::abs(liquid_phi_plus.data[i]) >
        std::abs(liquid_phi_minus.data[i])) {
      liquid_phi.data[i] = liquid_phi_minus.data[i];
    } else {
      liquid_phi.data[i] = liquid_phi_plus.data[i];
    }
  }
}

// adjust particle radii
void Simulation::adjust_particle_radii() {
  for (auto &p : particles) {
    float local_phi_val = trilerp_scalar_field(liquid_phi, p.position);
    p.radius = glm::clamp(std::abs(local_phi_val), 0.1f * h, 0.5f * h);
  }
}

// compute gradient norm with godunov scheme
void Simulation::norm_gradient() {
  norm_grad.set(1);
  for (int i = 1; i < norm_grad.sx - 1; i++) {
    for (int j = 1; j < norm_grad.sy - 1; j++) {
      for (int k = 1; k < norm_grad.sz - 1; k++) {
        // value for choosing upwinding direction
        float a = sig(i, j, k);
        // dxn = dx negative = backwards upwinding
        // dxp = dx positive = forwards upwinding
        float dxn = (liquid_phi(i, j, k) - liquid_phi(i - 1, j, k)) / h;
        float dxp = (liquid_phi(i + 1, j, k) - liquid_phi(i, j, k)) / h;
        float dyn = (liquid_phi(i, j, k) - liquid_phi(i, j - 1, k)) / h;
        float dyp = (liquid_phi(i, j + 1, k) - liquid_phi(i, j, k)) / h;
        float dzn = (liquid_phi(i, j, k) - liquid_phi(i, j, k - 1)) / h;
        float dzp = (liquid_phi(i, j, k + 1) - liquid_phi(i, j, k)) / h;
        float phidx = 0;
        float phidy = 0;
        float phidz = 0;
        // phi_dx
        if (a >= 0) {
          dxn = (dxn > 0) ? dxn * dxn : 0;
          dxp = (dxp < 0) ? dxp * dxp : 0;
          phidx = std::max(dxn, dxp);
        } else {
          dxn = (dxn < 0) ? dxn * dxn : 0;
          dxp = (dxp > 0) ? dxp * dxp : 0;
          phidx = std::max(dxn, dxp);
        }
        // phi_dy
        if (a >= 0) {
          dyn = (dyn > 0) ? dyn * dyn : 0;
          dyp = (dyp < 0) ? dyp * dyp : 0;
          phidy = std::max(dyn, dyp);
        } else {
          dyn = (dyn < 0) ? dyn * dyn : 0;
          dyp = (dyp > 0) ? dyp * dyp : 0;
          phidy = std::max(dyn, dyp);
        }
        // phi_dz
        if (a >= 0) {
          dzn = (dzn > 0) ? dzn * dzn : 0;
          dzp = (dzp < 0) ? dzp * dzp : 0;
          phidz = std::max(dzp, dzn);
        } else {
          dzn = (dzn < 0) ? dzn * dzn : 0;
          dzp = (dzp > 0) ? dzp * dzp : 0;
          phidz = std::max(dzp, dzn);
        }
        norm_grad(i, j, k) = std::sqrt(phidx + phidy + phidz);
      }
    }
  }
}

void Simulation::reinitialize_phi() {
  // compute sigmoid function for phi_0
  sig = liquid_phi / (liquid_phi * liquid_phi +
                      h * h); // yes this component-wise math is allowed
  norm_gradient();

  float err = 0;
  float tol = 1e-2;
  int max_iter = 2000;
  float dt = 0.1f * h;
  for (int iter = 0; iter <= max_iter; iter++) {
    if (iter == max_iter)
      throw std::runtime_error("error: phi reinitialization did not converge");

    // compute updated phi values for one timestep
    liquid_phi = liquid_phi - ((sig * (norm_grad - 1)) * dt);
    // recompute gradient norms
    norm_gradient();
    // compute average error
    err = 0;
    for (int i = 0; i < norm_grad.size; i++) {
      err += std::abs(norm_grad.data[i] - 1.0f);
    }
    err /= static_cast<float>(norm_grad.size);
    // std::cout << "err:" << err << "\n";

    if (err < tol) {
      // std::cout << "phi succesfully reinitialized\n";
      break;
    }
  }
}

// particle advection with bounds checking
void Simulation::advect_particles(float dt) {
  for (auto &p : particles) {
    p.position = rk2(p.position, dt);
    // bounds checking
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, glm::vec3(0), index, coords);
    float solid_phi_val = solid_phi.trilerp(index, coords);
    if (solid_phi_val < 0) {
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
  for (int i = 0; i < liquid_phi.sx; i++) {
    for (int j = 0; j < liquid_phi.sy; j++) {
      for (int k = 0; k < liquid_phi.sz; k++) {
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
void Simulation::add_gravity(float dt) {
  if (not predefined_field)
    v -= 9.8 * dt;
}

void Simulation::enforce_boundaries() {
  for (int i = 0; i < solid_phi.sx; i++) {
    for (int j = 0; j < solid_phi.sy; j++) {
      for (int k = 0; k < solid_phi.sz; k++) {
        // if solid cell, set face velocities to 0
        if (solid_phi(i, j, k) <= 0) {
          liquid_phi(i, j, k) = std::max(liquid_phi(i, j, k), 0.5f * h);
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
void Simulation::project(float dt) {
  compute_divergence();
  solve_pressure(dt);
  apply_pressure_gradient(dt);
}

// compute negative divergence to be used in rhs of pressure solve
void Simulation::compute_divergence() {
  divergence.clear();
  float coef = -1.0f / h;
  for (int i = 0; i < divergence.sx; i++) {
    for (int j = 0; j < divergence.sy; j++) {
      for (int k = 0; k < divergence.sz; k++) {
        if (liquid_phi(i, j, k) <= 0)
          divergence(i, j, k) =
              coef * (u(i + 1, j, k) - u(i, j, k) + v(i, j + 1, k) -
                      v(i, j, k) + w(i, j, k + 1) - w(i, j, k));
      }
    }
  }
}

void Simulation::solve_pressure_helper(std::vector<Eigen::Triplet<double>> &tl,
                                       double &aii, float dt, int i, int j,
                                       int k, int i1, int j1, int k1) {
  if (solid_phi(i1, j1, k1) > 0) {
    aii += 1;
    if (liquid_phi(i1, j1, k1) <= 0) {
      double scale = dt / (density * h * h);
      tl.push_back(Eigen::Triplet<double>(fluid_index(i, j, k),
                                          fluid_index(i1, j1, k1), -scale));
    }
  }
}

void Simulation::solve_pressure(float dt) {
  int nf = 0; // total number of fluid cells

  // assign an index to each fluid cell
  fluid_index.clear();
  for (int i = 0; i < fluid_index.sx; i++) {
    for (int j = 0; j < fluid_index.sy; j++) {
      for (int k = 0; k < fluid_index.sz; k++) {
        if (liquid_phi(i, j, k) <= 0) {
          fluid_index(i, j, k) = nf;
          nf += 1;
        }
      }
    }
  }

  // populate triplets
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(nf * 7);
  for (int i = 0; i < divergence.sx; i++) {
    for (int j = 0; j < divergence.sy; j++) {
      for (int k = 0; k < divergence.sz; k++) {
        if (liquid_phi(i, j, k) <= 0) {
          int index = fluid_index(i, j, k);
          double aii = 0; // negative sum of nonsolid nbrs
          double scale = dt / (density * h * h);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i + 1, j, k);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i - 1, j, k);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j + 1, k);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j - 1, k);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j, k + 1);
          solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j, k - 1);
          tripletList.push_back(
              Eigen::Triplet<double>(index, index, aii * scale));
        }
      }
    }
  }
  // construct A from triplets
  Eigen::SparseMatrix<double> A(nf, nf);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  // construct b
  Eigen::VectorXd b(nf);
  for (int i = 0; i < divergence.sx; i++) {
    for (int j = 0; j < divergence.sy; j++) {
      for (int k = 0; k < divergence.sz; k++) {
        if (liquid_phi(i, j, k) <= 0) {
          b(fluid_index(i, j, k)) = divergence(i, j, k);
        }
      }
    }
  }

  // solve Ax=b
  Eigen::VectorXd x(nf);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  x = solver.solve(b);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("pressure not solved");
  }

  // set pressure
  pressure.clear();
  for (int i = 0; i < fluid_index.sx; i++) {
    for (int j = 0; j < fluid_index.sy; j++) {
      for (int k = 0; k < fluid_index.sz; k++) {
        if (liquid_phi(i, j, k) <= 0) {
          pressure(i, j, k) = x(fluid_index(i, j, k));
        }
      }
    }
  }
}

// apply pressure update with ghost fluid method
// FIXME add variable densities
void Simulation::apply_pressure_gradient(float dt) {
  for (int i = 1; i < nx; i++) {
    for (int j = 1; j < ny; j++) {
      for (int k = 1; k < nz; k++) {
        // note: density is sampled at grid faces
        if (liquid_phi(i - 1, j, k) <= 0 || liquid_phi(i, j, k) <= 0)
          u(i, j, k) -= (dt / (density * h)) *
                        (pressure(i, j, k) - pressure(i - 1, j, k));

        if (liquid_phi(i, j - 1, k) <= 0 || liquid_phi(i, j, k) <= 0)
          v(i, j, k) -= (dt / (density * h)) *
                        (pressure(i, j, k) - pressure(i, j - 1, k));

        if (liquid_phi(i, j, k - 1) <= 0 || liquid_phi(i, j, k) <= 0)
          w(i, j, k) -= (dt / (density * h)) *
                        (pressure(i, j, k) - pressure(i, j, k - 1));
      }
    }
  }
}

// advance the simulation by a given time
// the time will be split up into substeps
// to satisfy the CFL condition
void Simulation::advance(float dt) {
  float t = 0.0f;
  float substep;
  while (t < dt) {
    substep = 0.5f * CFL();
    if (t + substep >= dt) {
      substep = dt - t;
    }
    t += substep;
    // particle level set method
    reseed_particles();
    advect_particles(substep);
    advect_phi(substep);
    correct_levelset();
    reinitialize_phi();
    correct_levelset();
    adjust_particle_radii();

    advect_velocity(substep);
    add_gravity(substep);
    enforce_boundaries();
    project(substep);
    enforce_boundaries();

    // extrapolate();
  }
}