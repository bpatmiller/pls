#pragma once

#include "array3.h"

#include <eigen3/Eigen/SparseCore>

#define AIR_CELL 0
#define FLUID_CELL 1
#define SOLID_CELL 2

class Grid {
public:
  float gravity = -9.8f; // force of gravity

  float lx, ly, lz; // total length of grid
  float nx, ny, nz; // number of cells per dimension
  float h;          // size of each cell
  // multiphase components
  std::vector<float> densities;   // per fluid densities
  std::vector<float> viscosities; // per fluid viscosities
  Array3<glm::vec4> phis; // signed distances from fluid, split into components
  Array3<glm::vec4> markers; // per fluid marker
  Array3f phi;               // minimum of each fluid phi
  Array3i ids;               // id of fluid at each grid cell (based on min phi)
  int fluid_count = 0;

  Array3f u, v, w;       // velocities sampled at cell faces
  Array3f u_w, v_w, w_w; // keep track of the total weight of particles
  Array3f du, dv, dw;    // saved velocities for flip
  Array3i marker;        // designates air, fluid, solid
  Array3d pressure;      // pressure at each grid cell
  Array3d r;             // divergence at each grid cell
  Array3i fl_index;      // gives each fluid cell an index (used for poission
                         // construction)
  Array3i pc;            // count used for reseeding

  // used for pressure solve
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd x, b;

  Grid() {}

  // adds a fluid to the grid, returns its id
  int add_fluid(float d, float v) {
    fluid_count++;
    densities.push_back(d);
    viscosities.push_back(v);
    return fluid_count - 1;
  }

  void init(float lx_, int nx_, int ny_, int nz_) {
    // number of cells
    nx = nx_;
    ny = ny_;
    nz = nz_;
    // cell size
    h = lx_ / static_cast<float>(nx);
    // grid side lengths
    lx = lx_;
    ly = h * static_cast<float>(ny);
    lz = h * static_cast<float>(nz);

    // velocities
    u.init(nx + 1, ny, nz);
    v.init(nx, ny + 1, nz);
    w.init(nx, ny, nz + 1);
    u_w.init(nx + 1, ny, nz);
    v_w.init(nx, ny + 1, nz);
    w_w.init(nx, ny, nz + 1);
    du.init(nx + 1, ny, nz);
    dv.init(nx, ny + 1, nz);
    dw.init(nx, ny, nz + 1);

    marker.init(nx, ny, nz);
    phis.init(nx, ny, nz);
    markers.init(nx, ny, nz);
    phi.init(nx, ny, nz);
    ids.init(nx, ny, nz);
    pressure.init(nx, ny, nz);
    r.init(nx, ny, nz);
    fl_index.init(nx, ny, nz);
    pc.init(nx, ny, nz);
  }

  void reset() {
    fluid_count = 0;
    u.clear();
    v.clear();
    w.clear();
    du.clear();
    dv.clear();
    dw.clear();
    u_w.clear();
    v_w.clear();
    w_w.clear();

    for (int i = 0; i < marker.size; i++) {
      if (marker.data[i] != SOLID_CELL)
        marker.data[i] = 0;
    }
    phis.clear();
    markers.clear();
    phi.clear();
    ids.clear();
    pressure.clear();
    r.clear();
    fl_index.clear();
    pc.clear();
  }

  float CFL();
  void add_gravity(float dt);
  void sweep_phi(int id);
  void solve_phi(float p, float q, float r, float &c);
  void extend_velocity();
  void sweep_velocity();
  void sweep_velocity_boundary(Array3f &arr);
  void sweep_u(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_v(int i0, int i1, int j0, int j1, int k0, int k1);
  void sweep_w(int i0, int i1, int j0, int j1, int k0, int k1);
  void enforce_boundary();

  void save_velocity();

  void project(float dt);
  void compute_divergence();
  void form_poisson(float dt);
  void apply_poisson(Array3d &x, Array3d &y);
  void form_preconditioner();
  void apply_preconditioner(Array3d &x, Array3d &y, Array3d &m);
  void solve_pressure();
  void add_pressure_gradient(float dt);

  void solve_x_helper(std::vector<Eigen::Triplet<double>> &tl, double &aii,
                      float dt, int i, int j, int k, int i1, int j1, int k1);
  void solve_x(float dt);
};