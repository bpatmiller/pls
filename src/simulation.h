#pragma once

#include "Particle.h"
#include "marchingcubes.h"
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
  // phi data
  Array3f solid_phi, liquid_phi; // solid and fluid signed distances
  int fluid_cell_count;          // for visualization purposes
  // velocity data
  Array3f u, v, w;    // velocity field sampled at grid faces
  Array3f tu, tv, tw; // for storing prior velocity field
  Array3f u_weight, v_weight, w_weight;
  Array3f pressure;
  // particle data
  float particle_radius;
  std::vector<Particle> particles;

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
    particle_radius = (float)(h * 1.01 * std::sqrt(3.0f) * 0.5f);
    // initialize arrays
    solid_phi.init(nx, ny, nz);
    liquid_phi.init(nx, ny, nz);
    pressure.init(nx, ny, nz);

    u.init(nx + 1, ny, nz);
    v.init(nx, ny + 1, nz);
    w.init(nx, ny, nz + 1);
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
  // auxillary methods
  void intialize_boundaries();
  void step_frame(float time);
  void advance(float dt);
  glm::vec3 rk2(glm::vec3 position, float dt);
  // simulation methods
  void advect_particles(float dt);
  void advect(float dt);
  void add_gravity(float dt);
  void project(float dt);
  // particle functions
  void position_to_grid(glm::vec3 p, glm::vec3 offset, glm::ivec3 &index,
                        glm::vec3 &coords);
  glm::vec3 trilerp_uvw(glm::vec3 p);

  void compute_phi();
};