#include "simulation.h"

// ----------------
// simulation steps
// ----------------

// returns a safe timestep
float Simulation::CFL() {
  float max_v =
      std::sqrt(std::pow(u.infnorm(), 2.0f) + std::pow(w.infnorm(), 2.0f) +
                std::pow(v.infnorm(), 2.0f));
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
bool Simulation::check_corners(Array3f &phi, int i, int j, int k) {
  for (int x = 0; x <= 1; x++) {
    for (int y = 0; y <= 1; y++) {
      for (int z = 0; z <= 1; z++) {
        float phi_val = phi(i + x, j + y, k + z);
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
  for (auto &f : fluids) {
    initialize_particles(f);
  }
}

void Simulation::initialize_particles(Fluid &fluid) {
  // count how many particles are in each cell
  fluid.particle_count.clear();
  for (auto &p : fluid.particles) {
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, glm::vec3(0), index, coords);
    fluid.particle_count(index) += 1;
  }

  // "top off" remaining cells with small phi values
  for (int i = 1; i < fluid.phi.sx - 2; i++) {
    for (int j = 1; j < fluid.phi.sy - 2; j++) {
      for (int k = 1; k < fluid.phi.sz - 2; k++) {
        if (check_corners(fluid.phi, i, j, k)) {
          glm::vec3 base_position((i + 0.5f) * h, (j + 0.5f) * h,
                                  (k + 0.5f) * h);
          for (int count = fluid.particle_count(i, j, k); count < 64; count++) {
            glm::vec3 position =
                base_position + glm::linearRand(glm::vec3(0.0f), glm::vec3(h));
            // attraction step
            float initial_phi_val = trilerp_scalar_field(fluid.phi, position);
            float phi_goal =
                (initial_phi_val > 0)
                    ? glm::clamp(initial_phi_val, 0.1f * h, 1.0f * h)
                    : glm::clamp(initial_phi_val, -1.0f * h, -0.1f * h);
            glm::vec3 normal = glm::normalize(fluid.phi.gradlerp(
                glm::ivec3(i, j, k), (position - base_position) / h, h));
            glm::vec3 new_position =
                glm::clamp(position + (phi_goal - initial_phi_val) * normal,
                           1.0001f * h, (nx - 1.0001f) * h);
            float new_phi_val = trilerp_scalar_field(fluid.phi, new_position);

            // add particle with position, starting phi, and radius
            Particle p(new_position, new_phi_val,
                       glm::clamp(std::abs(new_phi_val), 0.1f * h, 0.5f * h));
            fluid.particles.push_back(p);
            fluid.particle_count(i, j, k) += 1;
          }
        }
      }
    }
  }
}

void Simulation::remove_particles() {
  for (auto &f : fluids) {
    remove_particles(f);
  }
}

// remove all extraneous particles
void Simulation::remove_particles(Fluid &fluid) {
  // mark all distant particles to be deleted
  for (auto &p : fluid.particles) {
    float phi_val = trilerp_scalar_field(fluid.phi, p.position);
    p.valid = (std::abs(phi_val) < 2.0f * h);
  }
  // remove from particle vector
  fluid.particles.erase(
      std::remove_if(fluid.particles.begin(), fluid.particles.end(),
                     [](Particle const &p) { return !p.valid; }),
      fluid.particles.end());
}

void Simulation::correct_levelset() {
  for (auto f : fluids) {
    correct_levelset(f);
  }
}

void Simulation::correct_levelset(Fluid &fluid) {
  fluid.phi_minus = fluid.phi;
  fluid.phi_plus = fluid.phi;

  for (auto &p : fluid.particles) {
    // determine if the particle is escaped
    float local_phi_val = trilerp_scalar_field(fluid.phi, p.position);
    if (p.starting_phi * local_phi_val < 0 &&
        std::abs(local_phi_val) > p.radius) {
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
              fluid.phi_plus(index.x + i, index.y + j, index.z + k) = std::max(
                  phi_p, fluid.phi_plus(index.x + i, index.y + j, index.z + k));
            }
            // rebuild phi_-
            else {
              fluid.phi_minus(index.x + i, index.y + j, index.z + k) =
                  std::min(phi_p, fluid.phi_minus(index.x + i, index.y + j,
                                                  index.z + k));
            }
          }
        }
      }
    }
  }

  // merge phi_+ and phi_-
  // equation 13 in PLS
  for (int i = 0; i < fluid.phi.size; i++) {
    if (std::abs(fluid.phi_plus.data[i]) > std::abs(fluid.phi_minus.data[i])) {
      fluid.phi.data[i] = fluid.phi_minus.data[i];
    } else {
      fluid.phi.data[i] = fluid.phi_plus.data[i];
    }
  }
}

void Simulation::adjust_particle_radii() {
  for (auto &f : fluids) {
    adjust_particle_radii(f);
  }
}

// adjust particle radii
void Simulation::adjust_particle_radii(Fluid &fluid) {
  for (auto &p : fluid.particles) {
    float local_phi_val = trilerp_scalar_field(fluid.phi, p.position);
    p.radius = glm::clamp(std::abs(local_phi_val), 0.1f * h, 0.5f * h);
  }
}

// compute gradient norm with godunov scheme
// FIXME reimpliment this in a much cleaner way
// specifically ignore edge cases and use the max(min max)
// formulation
void Simulation::norm_gradient(Fluid &fluid) {
  float limit = (nx + ny + nz) * h;
  fluid.phi.clamp(-limit, limit);
  fluid.norm_grad.set(1.0f);

  // non-boundary cells
  for (int i = 1; i < nx - 1; i++) {
    for (int j = 1; j < ny - 1; j++) {
      for (int k = 1; k < nz - 1; k++) {
        // value for choosing upwinding direction
        float a = fluid.sig(i, j, k);
        // dxn = dx negative = backwards upwinding
        // dxp = dx positive = forwards upwinding

        float phidx = 0;
        float phidy = 0;
        float phidz = 0;
        // phi_dx

        float dxn = (fluid.phi(i, j, k) - fluid.phi(i - 1, j, k)) / h;
        float dxp = (fluid.phi(i + 1, j, k) - fluid.phi(i, j, k)) / h;
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

        float dyn = (fluid.phi(i, j, k) - fluid.phi(i, j - 1, k)) / h;
        float dyp = (fluid.phi(i, j + 1, k) - fluid.phi(i, j, k)) / h;
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
        float dzn = (fluid.phi(i, j, k) - fluid.phi(i, j, k - 1)) / h;
        float dzp = (fluid.phi(i, j, k + 1) - fluid.phi(i, j, k)) / h;

        if (a >= 0) {
          dzn = (dzn > 0) ? dzn * dzn : 0;
          dzp = (dzp < 0) ? dzp * dzp : 0;
          phidz = std::max(dzp, dzn);
        } else {
          dzn = (dzn < 0) ? dzn * dzn : 0;
          dzp = (dzp > 0) ? dzp * dzp : 0;
          phidz = std::max(dzp, dzn);
        }

        fluid.norm_grad(i, j, k) = std::sqrt(phidx + phidy + phidz);
      }
    }
  }
}

void Simulation::reinitialize_phi() {
  for (auto &f : fluids) {
    reinitialize_phi(f);
  }
}

void Simulation::reinitialize_phi(Fluid &fluid) {
  // compute sigmoid function for phi_0
  for (int i = 0; i < fluid.sig.size; i++) {
    fluid.sig.data[i] = fluid.phi.data[i] /
                        std::sqrt(std::pow(fluid.phi.data[i], 2.0f) + h * h);
  }
  norm_gradient(fluid);

  float err = 0;
  float tol = 1e-1f;
  int max_iter = 50;
  float dt = 0.05f * h;
  for (int iter = 0; iter <= max_iter; iter++) {
    if (iter == max_iter)
      // std::cout << "error: phi reinitialization did not converge\n";
      std::cout << "error: phi reinit did not converge\n";

    // compute updated phi values for one timestep
    // liquid_phi = liquid_phi - ((sig * (norm_grad - 1.0f)) * dt);
    for (int i = 1; i < nx - 1; i++) {
      for (int j = 1; j < ny - 1; j++) {
        for (int k = 1; k < nz - 1; k++) {
          fluid.phi(i, j, k) -=
              fluid.sig(i, j, k) * (fluid.norm_grad(i, j, k) - 1.0f) * dt;
        }
      }
    }
    // recompute gradient norms
    norm_gradient(fluid);
    // compute average error
    err = 0;
    for (int i = 0; i < fluid.norm_grad.size; i++) {
      err += std::abs(fluid.norm_grad.data[i] - 1.0f);
    }
    err /= static_cast<float>(fluid.norm_grad.size);
    // std::cout << "err:" << err << "\n";

    if (err < tol) {
      // std::printf(
      //     "phi successfully reinitialized with error %f < %f in %i
      //     iterations\n", err, tol, iter);
      break;
    }
  }
}

void Simulation::advect_particles(float dt) {
  for (auto &f : fluids) {
    advect_particles(f, dt);
  }
}

// particle advection with bounds checking
void Simulation::advect_particles(Fluid &fluid, float dt) {
  for (auto &p : fluid.particles) {
    p.position = rk2(p.position, dt);
    // data for bounds checking
    glm::vec3 coords;
    glm::ivec3 index;
    position_to_grid(p.position, glm::vec3(0), index, coords);
    float solid_phi_val = solid_phi.trilerp(index, coords);
    // only if the particle is a NEGATIVE marker particle
    if (solid_phi_val < 0 && p.starting_phi < 0) {
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
  float min_x = h * 1.001f;
  float max_x = h * ((float)nx - 1.001f);

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

void Simulation::advect_phi(float dt) {
  for (auto &f : fluids) {
    advect_phi(f, dt);
  }
}

// advect phi with equation phi_t = - V dot grad phi
void Simulation::advect_phi(Fluid &fluid, float dt) {
  fluid.phi_copy.clear();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        glm::vec3 position((i + 0.5f) * h, (j + 0.5f) * h, (k + 0.5f) * h);

        glm::vec3 velocity = trilerp_uvw(position);
        glm::vec3 phi_grad = fluid.phi.upwind_gradient(i, j, k, h, velocity);

        fluid.phi_copy(i, j, k) =
            fluid.phi(i, j, k) - dt * glm::dot(velocity, phi_grad);
      }
    }
  }
  fluid.phi_copy.copy_to(fluid.phi);
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
void Simulation::add_gravity(float dt) { v -= 9.8f * dt; }

void Simulation::enforce_boundaries() {
  for (int i = 0; i < solid_phi.sx; i++) {
    for (int j = 0; j < solid_phi.sy; j++) {
      for (int k = 0; k < solid_phi.sz; k++) {
        // if solid cell, set face velocities to 0
        if (solid_phi(i, j, k) <= 0) {
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

  // make sure we dont have fluid overlapping with solid
  for (auto &f : fluids) {
    for (int i = 0; i < solid_phi.size; i++) {
      if (solid_phi.data[i] <= 0)
        f.phi.data[i] = std::max(f.phi.data[i], 0.5f * h);
    }
  }
}

void Simulation::project_phi() {
  for (i = 0; i < (nx * ny * nz); i++) {
    // get the lowest two phi values, and subtract the avg
    // from all phi values
    // NOTE: assumes at least 2 phase flow
    if (solid_phi.data[i] > 0) {
      // find the lowest 2
      float min1 = nx * ny * nz;
      float min2 = nx * ny * nz;
      for (auto & f : fluids) {
        if (f.phi.data[i] < min1) {
          min2 = min1;
          min1 = f.phi.data[i];
        } else if (f.phi.data[i] < min2) {
          min2 = f.phi.data[i];
        }
        // now subtract the average
        float avg = (min1 + min2) * 0.5f;
        for (auto & f : fluids) {
          f.phi.data[i] -= avg;
        }
      }
    }
  }
}

// project the velocity field onto its divergence-free part
void Simulation::project(float dt) {
  compute_divergence();
  // solve_pressure(dt);
  // apply_pressure_gradient(dt);
}

// compute negative divergence to be used in rhs of pressure solve
void Simulation::compute_divergence() {
  divergence.clear();
  for (int i = 0; i < divergence.sx; i++) {
    for (int j = 0; j < divergence.sy; j++) {
      for (int k = 0; k < divergence.sz; k++) {
        divergence(i, j, k) = u(i + 1, j, k) - u(i, j, k) + v(i, j + 1, k) -
                              v(i, j, k) + w(i, j, k + 1) - w(i, j, k);
      }
    }
  }
}

// void Simulation::solve_pressure_helper(std::vector<Eigen::Triplet<double>>
// &tl,
//                                        double &aii, float dt, int i, int j,
//                                        int k, int i1, int j1, int k1) {
//   if (solid_phi(i1, j1, k1) > 0) {
//     aii += 1;
//     if (liquid_phi(i1, j1, k1) <= 0) {
//       double scale = dt / (density_1 * h * h);
//       tl.push_back(Eigen::Triplet<double>(fluid_index(i, j, k),
//                                           fluid_index(i1, j1, k1), -scale));
//     }
//   }
// }

// void Simulation::solve_pressure(float dt) {
//   int nf = 0; // total number of fluid cells

//   // assign an index to each fluid cell
//   fluid_index.clear();
//   for (int i = 0; i < fluid_index.sx; i++) {
//     for (int j = 0; j < fluid_index.sy; j++) {
//       for (int k = 0; k < fluid_index.sz; k++) {
//         if (liquid_phi(i, j, k) <= 0) {
//           fluid_index(i, j, k) = nf;
//           nf += 1;
//         }
//       }
//     }
//   }

//   // populate triplets
//   std::vector<Eigen::Triplet<double>> tripletList;
//   tripletList.reserve(nf * 7);
//   for (int i = 0; i < divergence.sx; i++) {
//     for (int j = 0; j < divergence.sy; j++) {
//       for (int k = 0; k < divergence.sz; k++) {
//         if (liquid_phi(i, j, k) <= 0) {
//           int index = fluid_index(i, j, k);
//           double aii = 0; // negative sum of nonsolid nbrs
//           double scale = dt / (density_1 * h * h);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i + 1, j, k);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i - 1, j, k);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j + 1, k);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j - 1, k);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j, k + 1);
//           solve_pressure_helper(tripletList, aii, dt, i, j, k, i, j, k - 1);
//           tripletList.push_back(
//               Eigen::Triplet<double>(index, index, aii * scale));
//         }
//       }
//     }
//   }
//   // construct A from triplets
//   Eigen::SparseMatrix<double> A(nf, nf);
//   A.setFromTriplets(tripletList.begin(), tripletList.end());

//   // construct b
//   Eigen::VectorXd b(nf);
//   for (int i = 0; i < divergence.sx; i++) {
//     for (int j = 0; j < divergence.sy; j++) {
//       for (int k = 0; k < divergence.sz; k++) {
//         if (liquid_phi(i, j, k) <= 0) {
//           b(fluid_index(i, j, k)) = (-1.0 / h) * divergence(i, j, k);
//         }
//       }
//     }
//   }

//   // solve Ax=b
//   Eigen::VectorXd x(nf);
//   Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
//   solver.compute(A);
//   x = solver.solve(b);

//   if (solver.info() != Eigen::Success) {
//     throw std::runtime_error("pressure not solved");
//   }

//   // set pressure
//   pressure.clear();
//   for (int i = 0; i < fluid_index.sx; i++) {
//     for (int j = 0; j < fluid_index.sy; j++) {
//       for (int k = 0; k < fluid_index.sz; k++) {
//         if (liquid_phi(i, j, k) <= 0) {
//           pressure(i, j, k) = x(fluid_index(i, j, k));
//         }
//       }
//     }
//   }
// }

// // apply pressure update
// // TODO use ghost fluid method
// // FIXME add variable densities
// void Simulation::apply_pressure_gradient(float dt) {
//   for (int i = 1; i < nx; i++) {
//     for (int j = 1; j < ny; j++) {
//       for (int k = 1; k < nz; k++) {
//         // note: density is sampled at grid faces
//         if ((solid_phi(i, j, k) > 0 && solid_phi(i - 1, j, k) > 0) &&
//             (liquid_phi(i - 1, j, k) <= 0 || liquid_phi(i, j, k) <= 0)) {
//           u(i, j, k) -= (dt / (density_1 * h)) *
//                         (pressure(i, j, k) - pressure(i - 1, j, k));
//         }

//         if ((solid_phi(i, j, k) > 0 && solid_phi(i, j - 1, k) > 0) &&
//             (liquid_phi(i, j - 1, k) <= 0 || liquid_phi(i, j, k) <= 0)) {
//           v(i, j, k) -= (dt / (density_1 * h)) *
//                         (pressure(i, j, k) - pressure(i, j - 1, k));
//         }

//         if ((solid_phi(i, j, k) > 0 && solid_phi(i, j, k - 1) > 0) &&
//             (liquid_phi(i, j, k - 1) <= 0 || liquid_phi(i, j, k) <= 0)) {
//           w(i, j, k) -= (dt / (density_1 * h)) *
//                         (pressure(i, j, k) - pressure(i, j, k - 1));
//         }
//       }
//     }
//   }
// }

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

    // particle level set method
    advect_phi(substep);
    advect_particles(substep);
    correct_levelset();
    reinitialize_phi();
    correct_levelset();
    adjust_particle_radii();
    reseed_particles();

    project_phi();

    advect_velocity(substep);
    add_gravity(substep);

    enforce_boundaries();
    project(substep); // TODO
  }
}