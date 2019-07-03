#pragma once

#include "Particle.h"
#include "marchingcubes.h"
#include "util.h"
#include <algorithm>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseCore>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>
#include <glm/gtx/common.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <vector>

#define U_OFFSET glm::vec3(0.5f, 0, 0)
#define V_OFFSET glm::vec3(0, 0.5f, 0)
#define W_OFFSET glm::vec3(0, 0, 0.5f)
#define CENTER_OFFSET glm::vec3(0.5f, 0.5f, 0.5f)

class Simulation {
public:
  int nx, ny, nz;   // grid resolution
  float lx, ly, lz; // container width
  float h;          // cell width
  float density_1 = 1.0f;
  float density_2 = 2.0f;
  // phi data
  Array3f solid_phi, liquid_phi, phi_copy; // solid and fluid signed distances
  Array3f sig, norm_grad; // sigmoid smoothing function and norm of gradient
  Array3f liquid_phi_plus, liquid_phi_minus;
  // velocity data
  Array3f u, v, w;             // velocity field sampled at grid faces
  Array3f tu, tv, tw;          // for storing prior velocity field
  Array3f u_vol, v_vol, w_vol; // volume fractions for pressure solve
  Array3f u_weight, v_weight, w_weight;
  bool predefined_field = false;
  // pressure data
  Array3d pressure;
  Array3d divergence;
  Array3i fluid_index;
  // particle data
  std::vector<Particle> particles;
  Array3i particle_count;
  int reseed_counter;

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x, b;

  Simulation(){};

  void init(float lx_, int nx_, int ny_, int nz_) {
    // set up dimension variables
    nx = nx_;
    ny = ny_;
    nz = nz_;
    lx = lx_;
    h = lx / static_cast<float>(nx);
    ly = h * ny;
    lz = h * nz;
    reseed_counter = 0;
    // initialize arrays
    solid_phi.init(nx, ny, nz);
    liquid_phi.init(nx, ny, nz);
    liquid_phi_minus.init(nx, ny, nz);
    liquid_phi_plus.init(nx, ny, nz);
    phi_copy.init(nx, ny, nz);
    sig.init(nx, ny, nz);
    norm_grad.init(nx, ny, nz);
    pressure.init(nx, ny, nz);
    divergence.init(nx, ny, nz);
    fluid_index.init(nx, ny, nz);

    u.init(nx + 1, ny, nz);
    u_vol.init(nx + 1, ny, nz);
    v.init(nx, ny + 1, nz);
    v_vol.init(nx, ny + 1, nz);
    w.init(nx, ny, nz + 1);
    w_vol.init(nx, ny, nx + 1);
    tu.init(nx + 1, ny, nz);
    tv.init(nx, ny + 1, nz);
    tw.init(nx, ny, nz + 1);
    u_weight.init(nx + 1, ny, nz);
    v_weight.init(nx, ny + 1, nz);
    w_weight.init(nx, ny, nz + 1);

    particle_count.init(nx, ny, nz);

    // set up solid boundaries
    intialize_boundaries();
  }

  float CFL();
  // initialization methods
  void add_sphere_phi();
  void add_cylinder_phi();
  void set_vortex_field();
  // auxillary methods
  void intialize_boundaries();
  void step_frame(float time);
  void advance(float dt);
  glm::vec3 rk2(glm::vec3 position, float dt);
  float trilerp_scalar_field(Array3f &field, glm::vec3 position);
  // simulation methods
  void advect_phi(float dt);
  void advect_velocity(float dt);
  void reinitialize_phi();
  void add_gravity(float dt);
  void enforce_boundaries();
  void project(float dt);
  void compute_divergence();
  void solve_pressure(float dt);
  void solve_pressure_helper(std::vector<Eigen::Triplet<double>> &tl,
                             double &aii, float dt, int i, int j, int k, int i1,
                             int j1, int k1);
  void apply_pressure_gradient(float dt);
  // particle functions
  void position_to_grid(glm::vec3 p, glm::vec3 offset, glm::ivec3 &index,
                        glm::vec3 &coords);
  glm::vec3 trilerp_uvw(glm::vec3 p);
  glm::vec3 trilerp_tutvtw(glm::vec3 p);
  bool check_corners(int i, int j, int k);
  void reseed_particles();
  void initialize_particles();
  void remove_particles();
  void advect_particles(float dt);
  void correct_levelset();
  void adjust_particle_radii();
  void compute_volume_fractions();
  void compute_volume_fractions_arr(Array3f &vol, Array3f &field,
                                    glm::vec3 offset);

  // auxillary level set methods
  void norm_gradient();
  // velocity extension
  void extend_velocity();
  void sweep_velocity();
  void sweep_velocity_boundary(Array3f &arr);
  void sweep_u(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_v(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_w(int i0, int i1, int j0, int j1, int k0, int k1);
};