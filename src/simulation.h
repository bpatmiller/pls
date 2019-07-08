#pragma once

#include "Particle.h"
#include "geom.h"
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

class Fluid {
public:
  float density;
  Array3f phi, phi_copy, phi_plus, phi_minus;
  Array3f sig, norm_grad;
  Array3i particle_count;
  std::vector<Particle> particles;

  Fluid(float density_, int nx, int ny, int nz) {
    density = density_;
    phi.init(nx, ny, nz);
    phi_copy.init(nx, ny, nz);
    phi_plus.init(nx, ny, nz);
    phi_minus.init(nx, ny, nz);
    sig.init(nx, ny, nz);
    norm_grad.init(nx, ny, nz);
    particle_count.init(nx, ny, nz);
  }
};

class Simulation {
public:
  int nx, ny, nz;   // grid resolution
  float lx, ly, lz; // container width
  float h;          // cell width
  // solid phi data
  Array3f solid_phi;
  Array3f liquid_phi;
  // velocity data
  Array3f u, v, w;             // velocity field sampled at grid faces
  Array3f tu, tv, tw;          // for storing prior velocity field
  Array3f u_vol, v_vol, w_vol; // volume fractions for pressure solve
  Array3f u_weight, v_weight, w_weight;
  bool predefined_field = false;
  // pressure data
  Array3d pressure;
  Array3d divergence;
  int reseed_counter;
  Array3i grid_ids;
  Array3i fl_index;
  float fixed_density = 1.0f;

  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x, b;

  std::vector<Fluid> fluids;

  Simulation(){};

  void init();

  Simulation(float lx_, int nx_, int ny_, int nz_) {
    init(lx_, nx_, ny_, nz_);
  };

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
    pressure.init(nx, ny, nz);
    divergence.init(nx, ny, nz);
    grid_ids.init(nx, ny, nz);
    fl_index.init(nx, ny, nz);

    u.init(nx + 1, ny, nz);
    u_vol.init(nx + 1, ny, nz);
    v.init(nx, ny + 1, nz);
    v_vol.init(nx, ny + 1, nz);
    w.init(nx, ny, nz + 1);
    w_vol.init(nx, ny, nz + 1);
    tu.init(nx + 1, ny, nz);
    tv.init(nx, ny + 1, nz);
    tw.init(nx, ny, nz + 1);
    u_weight.init(nx + 1, ny, nz);
    v_weight.init(nx, ny + 1, nz);
    w_weight.init(nx, ny, nz + 1);

    // set up solid boundaries
    intialize_boundaries();
  }

  float CFL();
  // initialization methods
  void add_sphere_phi();
  void add_two_spheres();
  void add_two_walls();
  void add_pool_phi();
  void add_pool_drop();
  void add_cylinder_phi();
  void set_vortex_field();
  void add_fluid(float density_);
  // auxillary methods
  void intialize_boundaries();
  void step_frame(float time);
  void advance(float dt);
  glm::vec3 rk2(glm::vec3 position, float dt);
  float trilerp_scalar_field(Array3f &field, glm::vec3 position);
  int ijk_to_index(int i, int j, int k);
  void norm_gradient(Fluid &fluid);
  // particle functions
  void position_to_grid(glm::vec3 p, glm::vec3 offset, glm::ivec3 &index,
                        glm::vec3 &coords);
  glm::vec3 trilerp_uvw(glm::vec3 p);
  glm::vec3 trilerp_tutvtw(glm::vec3 p);
  // multifluid-independent components
  void compute_divergence();
  void advect_velocity(float dt);
  void add_gravity(float dt);

  // Level set functions
  bool check_corners(Array3f &phi, int i, int j, int k);
  void reseed_particles();
  void initialize_particles();
  void initialize_particles(Fluid &fluid);
  void remove_particles();
  void remove_particles(Fluid &fluid);
  void correct_levelset();
  void correct_levelset(Fluid &fluid);
  void adjust_particle_radii();
  void adjust_particle_radii(Fluid &fluid);
  void reinitialize_phi();
  void reinitialize_phi(Fluid &fluid);
  void advect_particles(float dt);
  void advect_particles(Fluid &fluid, float dt);
  void advect_phi(float dt);
  void advect_phi(Fluid &fluid, float dt);

  void enforce_boundaries();
  void project_phi(); // from Multiple Interacting Liquids
  void project(float dt);

  // TODO
  void solve_pressure(float dt);
  void solve_pressure_helper(std::vector<Eigen::Triplet<double>> &tl,
                             double &aii, float dt, int i, int j, int k, int i1,
                             int j1, int k1);
  void apply_pressure_gradient(float dt);

  // velocity extrapolation
  void extend_velocity();
  void sweep_velocity();
  void sweep_velocity_boundary(Array3f &arr);
  void sweep_u(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_v(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_w(int i0, int i1, int j0, int j1, int k0, int k1);
  // io methods
  void export_mesh(float time);
  void generate_mesh(Fluid fluid, Mesh &mesh);
};