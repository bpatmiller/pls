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
        if (std::abs(phi_val) < 2.0f * h)
          return true;
      }
    }
  }
  return false;
}

void Simulation::reseed_particles() {
  if (reseed_counter % 5 == 0) {
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
  for (int i = 1; i < liquid_phi.sx - 2; i++) {
    for (int j = 1; j < liquid_phi.sy - 2; j++) {
      for (int k = 1; k < liquid_phi.sz - 2; k++) {
        if (check_corners(i, j, k)) {
          glm::vec3 base_position((i + 0.5f) * h, (j + 0.5f) * h,
                                  (k + 0.5f) * h);
          for (int count = particle_count(i, j, k); count < 64; count++) {
            glm::vec3 position =
                base_position + glm::linearRand(glm::vec3(0.0f), glm::vec3(h));
            // attraction step
            float initial_phi_val = trilerp_scalar_field(liquid_phi, position);
            float phi_goal =
                (initial_phi_val > 0)
                    ? glm::clamp(initial_phi_val, 0.1f * h, 2.0f * h)
                    : glm::clamp(initial_phi_val, -2.0f * h, -0.1f * h);
            glm::vec3 normal = glm::normalize(liquid_phi.gradlerp(
                glm::ivec3(i, j, k), (position - base_position) / h, h));
            glm::vec3 new_position =
                glm::clamp(position + (phi_goal - initial_phi_val) * normal,
                           1.0001f * h, (nx - 1.0001f) * h);
            float new_phi_val = trilerp_scalar_field(liquid_phi, new_position);

            // add particle with position, starting phi, and radius
            Particle p(new_position, new_phi_val,
                       glm::clamp(std::abs(new_phi_val), 0.1f * h, 0.5f * h));
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
  sig = liquid_phi / (liquid_phi * liquid_phi + h * h);
  norm_gradient();

  float err = 0;
  float tol = 1e-1f;
  int max_iter = 1000;
  float dt = 1e-1f * h;
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
  if (predefined_field)
    return;
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
    v -= 9.8f * dt;
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
  if (predefined_field)
    return;
  compute_divergence();
  solve_pressure(dt);
  apply_pressure_gradient(dt);
}

// compute negative divergence to be used in rhs of pressure solve
void Simulation::compute_divergence() {
  divergence.clear();
  for (int i = 0; i < divergence.sx; i++) {
    for (int j = 0; j < divergence.sy; j++) {
      for (int k = 0; k < divergence.sz; k++) {
        if (liquid_phi(i, j, k) <= 0)
          divergence(i, j, k) = u(i + 1, j, k) - u(i, j, k) + v(i, j + 1, k) -
                                v(i, j, k) + w(i, j, k + 1) - w(i, j, k);
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
      double scale = dt / (density_1 * h * h);
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
          double scale = dt / (density_1 * h * h);
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
          b(fluid_index(i, j, k)) = (-1.0 / h) * divergence(i, j, k);
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

// apply pressure update
// TODO use ghost fluid method
// FIXME add variable densities
void Simulation::apply_pressure_gradient(float dt) {
  for (int i = 2; i < nx; i++) {
    for (int j = 2; j < ny; j++) {
      for (int k = 2; k < nz; k++) {
        // note: density is sampled at grid faces
        if ((solid_phi(i, j, k) > 0 && solid_phi(i - 1, j, k) > 0) &&
            (liquid_phi(i - 1, j, k) <= 0 || liquid_phi(i, j, k) <= 0)) {
          u(i, j, k) -= (dt / (density_1 * h)) *
                        (pressure(i, j, k) - pressure(i - 1, j, k));
        }

        if ((solid_phi(i, j, k) > 0 && solid_phi(i, j - 1, k) > 0) &&
            (liquid_phi(i, j - 1, k) <= 0 || liquid_phi(i, j, k) <= 0)) {
          v(i, j, k) -= (dt / (density_1 * h)) *
                        (pressure(i, j, k) - pressure(i, j - 1, k));
        }

        if ((solid_phi(i, j, k) > 0 && solid_phi(i, j, k - 1) > 0) &&
            (liquid_phi(i, j, k - 1) <= 0 || liquid_phi(i, j, k) <= 0)) {
          w(i, j, k) -= (dt / (density_1 * h)) *
                        (pressure(i, j, k) - pressure(i, j, k - 1));
        }
      }
    }
  }
}

void Simulation::extend_velocity() {
  for (int i = 0; i < 8; i++) {
    sweep_velocity();
  }
}

void Simulation::sweep_velocity() {
  // U --------------------------------
  sweep_u(1, u.sx - 1, 1, u.sy - 1, 1, u.sz - 1);
  sweep_u(1, u.sx - 1, 1, u.sy - 1, u.sz - 2, 0);
  sweep_u(1, u.sx - 1, u.sy - 2, 0, 1, u.sz - 1);
  sweep_u(1, u.sx - 1, u.sy - 2, 0, u.sz - 2, 0);
  sweep_u(u.sx - 2, 0, 1, u.sy - 1, 1, u.sz - 1);
  sweep_u(u.sx - 2, 0, 1, u.sy - 1, u.sz - 2, 0);
  sweep_u(u.sx - 2, 0, u.sy - 2, 0, 1, u.sz - 1);
  sweep_u(u.sx - 2, 0, u.sy - 2, 0, u.sz - 2, 0);
  // set boundary cells
  sweep_velocity_boundary(u);

  // V --------------------------------
  sweep_v(1, v.sx - 1, 1, v.sy - 1, 1, v.sz - 1);
  sweep_v(1, v.sx - 1, 1, v.sy - 1, v.sz - 2, 0);
  sweep_v(1, v.sx - 1, v.sy - 2, 0, 1, v.sz - 1);
  sweep_v(1, v.sx - 1, v.sy - 2, 0, v.sz - 2, 0);
  sweep_v(v.sx - 2, 0, 1, v.sy - 1, 1, v.sz - 1);
  sweep_v(v.sx - 2, 0, 1, v.sy - 1, v.sz - 2, 0);
  sweep_v(v.sx - 2, 0, v.sy - 2, 0, 1, v.sz - 1);
  sweep_v(v.sx - 2, 0, v.sy - 2, 0, v.sz - 2, 0);
  // set boundary cells
  sweep_velocity_boundary(v);

  // W --------------------------------
  sweep_w(1, w.sx - 1, 1, w.sy - 1, 1, w.sz - 1);
  sweep_w(1, w.sx - 1, 1, w.sy - 1, w.sz - 2, 0);
  sweep_w(1, w.sx - 1, w.sy - 2, 0, 1, w.sz - 1);
  sweep_w(1, w.sx - 1, w.sy - 2, 0, w.sz - 2, 0);
  sweep_w(w.sx - 2, 0, 1, w.sy - 1, 1, w.sz - 1);
  sweep_w(w.sx - 2, 0, 1, w.sy - 1, w.sz - 2, 0);
  sweep_w(w.sx - 2, 0, w.sy - 2, 0, 1, w.sz - 1);
  sweep_w(w.sx - 2, 0, w.sy - 2, 0, w.sz - 2, 0);
  // set boundary cells
  sweep_velocity_boundary(w);
}

void Simulation::sweep_velocity_boundary(Array3f &arr) {
  // top and bottom
  for (int i = 0; i < arr.sx; i++) {
    for (int k = 0; k < arr.sz; k++) {
      arr(i, 0, k) = arr(i, 1, k);
      arr(i, arr.sy - 1, k) = arr(i, arr.sy - 2, k);
    }
  }
  // left and right
  for (int j = 0; j < arr.sy; j++) {
    for (int k = 0; k < arr.sz; k++) {
      arr(0, j, k) = arr(1, j, k);
      arr(arr.sx - 1, j, k) = arr(arr.sx - 2, j, k);
    }
  }
  // front and back
  for (int i = 0; i < arr.sx; i++) {
    for (int j = 0; j < arr.sy; j++) {
      arr(i, j, 0) = arr(i, j, 1);
      arr(i, j, arr.sz - 1) = arr(i, j, arr.sz - 2);
    }
  }
}

void Simulation::sweep_u(int i0, int i1, int j0, int j1, int k0, int k1) {
  int di = (i0 < i1) ? 1 : -1;
  int dj = (j0 < j1) ? 1 : -1;
  int dk = (k0 < k1) ? 1 : -1;

  float weight;

  for (int i = i0; i != i1; i += di) {
    for (int j = j0; j != j1; j += dj) {
      for (int k = k0; k != k1; k += dk) {
        if (liquid_phi(i - 1, j, k) > 0 && liquid_phi(i, j, k) > 0) {
          float dp = di * (liquid_phi(i, j, k) - liquid_phi(i - 1, j, k));
          if (dp < 0)
            continue;
          // avg y-dir phi change
          float dq =
              0.5 * (liquid_phi(i - 1, j, k) + liquid_phi(i, j, k) -
                     liquid_phi(i - 1, j - dj, k) - liquid_phi(i, j - dj, k));
          if (dq < 0)
            continue;
          // avg z-dir phi change
          float dr =
              0.5 * (liquid_phi(i - 1, j, k) + liquid_phi(i, j, k) -
                     liquid_phi(i - 1, j, k - dk) - liquid_phi(i, j, k - dk));
          if (dr < 0)
            continue;

          // weighted sum of other velocities
          if (dp + dq + dr == 0) {
            weight = 1.0f / 3.0f;
            u(i, j, k) =
                weight * (u(i - di, j, k) + u(i, j - dj, k) + u(i, j, k - dk));
          } else {
            weight = 1.0f / (dp + dq + dr);
            u(i, j, k) = dp * weight * u(i - di, j, k) +
                         dq * weight * u(i, j - dj, k) +
                         dr * weight * u(i, j, k - dk);
          }
        }
      }
    }
  }
}
void Simulation::sweep_v(int i0, int i1, int j0, int j1, int k0, int k1) {
  int di = (i0 < i1) ? 1 : -1;
  int dj = (j0 < j1) ? 1 : -1;
  int dk = (k0 < k1) ? 1 : -1;

  float weight;

  for (int i = i0; i != i1; i += di) {
    for (int j = j0; j != j1; j += dj) {
      for (int k = k0; k != k1; k += dk) {
        if (liquid_phi(i, j - 1, k) > 0 && liquid_phi(i, j, k) > 0) {
          float dq = dj * (liquid_phi(i, j, k) - liquid_phi(i, j - 1, k));
          if (dq < 0)
            continue;
          // avg x-dir phi change
          float dp =
              0.5 * (liquid_phi(i, j - 1, k) + liquid_phi(i, j, k) -
                     liquid_phi(i - di, j - 1, k) - liquid_phi(i - di, j, k));
          if (dp < 0)
            continue;
          // avg z-dir phi change
          float dr = 0.5 * (liquid_phi(i - 1, j - 1, k) + liquid_phi(i, j, k) -
                            liquid_phi(i - 1, j - 1, k - dk) -
                            liquid_phi(i, j, k - dk));
          if (dr < 0)
            continue;

          // weighted sum of other velocities
          if (dp + dq + dr == 0) {
            weight = 1.0f / 3.0f;
            v(i, j, k) =
                weight * (v(i - di, j, k) + v(i, j - dj, k) + v(i, j, k - dk));
          } else {
            weight = 1.0f / (dp + dq + dr);
            v(i, j, k) = dp * weight * v(i - di, j, k) +
                         dq * weight * v(i, j - dj, k) +
                         dr * weight * v(i, j, k - dk);
          }
        }
      }
    }
  }
}

void Simulation::sweep_w(int i0, int i1, int j0, int j1, int k0, int k1) {
  int di = (i0 < i1) ? 1 : -1;
  int dj = (j0 < j1) ? 1 : -1;
  int dk = (k0 < k1) ? 1 : -1;

  float weight;

  for (int i = i0; i != i1; i += di) {
    for (int j = j0; j != j1; j += dj) {
      for (int k = k0; k != k1; k += dk) {
        if (liquid_phi(i, j, k - 1) > 0 && liquid_phi(i, j, k) > 0) {
          float dr = dk * (liquid_phi(i, j, k) - liquid_phi(i, j, k - 1));
          if (dr < 0)
            continue;
          // avg y-dir phi change
          float dq = 0.5 * (liquid_phi(i, j - 1, k) + liquid_phi(i, j, k) -
                            liquid_phi(i, j - dj, k - dk) -
                            liquid_phi(i, j - 1, k - dk));
          if (dq < 0)
            continue;
          // avg x-dir phi change
          float dp = 0.5 * (liquid_phi(i - 1, j, k) + liquid_phi(i, j, k) -
                            liquid_phi(i - 1, j - 1, k - dk) -
                            liquid_phi(i - 1, j, k - dk));
          if (dp < 0)
            continue;

          // weighted sum of other velocities
          if (dp + dq + dr == 0) {
            weight = 1.0f / 3.0f;
            w(i, j, k) =
                weight * (w(i - di, j, k) + w(i, j - dj, k) + w(i, j, k - dk));
          } else {
            weight = 1.0f / (dp + dq + dr);
            w(i, j, k) = dp * weight * w(i - di, j, k) +
                         dq * weight * w(i, j - dj, k) +
                         dr * weight * w(i, j, k - dk);
          }
        }
      }
    }
  }
}

void Simulation::compute_volume_fractions_arr(Array3f &vol, Array3f &field,
                                              glm::vec3 offset) {
  float substep = h / 2.0f;
  float num_subcells = 8.0;
  // compute u-cell volume fractions
  for (int i = 1; i < field.sx - 1; i++) {
    for (int j = 1; j < field.sy - 1; j++) {
      for (int k = 1; k < field.sz - 1; k++) {
        float sum_fractions = 0;
        glm::vec3 base_position = glm::vec3(i * h, j * h, k * h) + offset;
        for (float subx = substep / 2.0f; subx < h; subx += substep) {
          for (float suby = substep / 2.0f; suby < h; suby += substep) {
            for (float subz = substep / 2.0f; subz < h; subz += substep) {
              glm::vec3 position = glm::vec3(subx, suby, subz) + base_position;
              float phi_val = trilerp_scalar_field(liquid_phi, position);
              // bridson 6.4
              sum_fractions += 0.5f - 0.5f * glm::clamp(4 * phi_val / h, -1.0f, 1.0f);
            }
          }
        }
        vol(i, j, k) = std::max(sum_fractions / num_subcells,1e-3f);
      }
    }
  }
}

void Simulation::compute_volume_fractions() {
  compute_volume_fractions_arr(u_vol, u, glm::vec3(-0.5f * h, 0.0f, 0.0f));
  compute_volume_fractions_arr(v_vol, v, glm::vec3(0.0f, -0.5f * h, 0.0f));
  compute_volume_fractions_arr(w_vol, w, glm::vec3(0.0f, 0.0f, -0.5f * h));
}

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
    advect_velocity(substep);
    add_gravity(substep);
    // particle level set method
    advect_phi(substep);
    advect_particles(substep);
    correct_levelset();
    reinitialize_phi();
    correct_levelset();
    adjust_particle_radii();
    reseed_particles();

    enforce_boundaries();
    project(substep);
    // enforce_boundaries();
    extend_velocity();
    // enforce_boundaries();
  }
}