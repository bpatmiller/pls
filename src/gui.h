#pragma once

#include <glad/glad.h>

#include "gl/program.h"
#include "gl/vao.h"
#include "simulation.h"
#include <GLFW/glfw3.h>
#include <glm/gtx/string_cast.hpp>
#include <map>

class GUI {
public:
  GUI(GLFWwindow *w) { window = w; }
  void init(float lx_, int nx_, int ny_, int nz_);
  void update(float t = -1, bool force = false);
  void update_camera();
  void create_sphere(float Radius, std::vector<glm::vec3> &s_vertices);

  // simulation data
  Simulation simulation;
  float timestep = 0.025f;
  float totaltime = 0.0f;
  // gl variables
  GLFWwindow *window;
  glm::ivec2 window_dims;

  // key controls
  std::map<int, bool> keyHeld;

  // mouse controls
  bool mouse_pressed = false;
  glm::vec2 mouse_pos = glm::vec2(-1, -1);
  glm::vec2 mouse_pos_prev = glm::vec2(-1, -1);
  glm::vec2 mouse_diff = glm::vec2(0, 0);

  // camera variables
  glm::quat orientation = glm::quat(glm::mat4(1));
  glm::vec3 eye = glm::vec3(0, 0, 4);
  glm::vec3 base_eye = glm::vec3(1, 1, 4);
  glm::vec3 focus = glm::vec3(1, 1, 1);
  glm::mat4 view_matrix;
  glm::mat4 projection_matrix;
  float pitch = 0;
  float yaw = 0;

  int levelset_to_draw = 0;
  bool draw_grid = true;
  bool display_phi = true;
  bool draw_particles = true;
  bool draw_velocity = false;
  bool dirty = true;
  // grid visualization
  Program grid_program;
  VAO grid_vao;
  std::vector<glm::uvec3> box_indices;
  std::vector<glm::vec4> grid_offsets;
  // marker particle visualization
  Program particle_program;
  VAO particle_vao;
  std::vector<glm::uvec3> sphere_indices;
  // velocity field visualization
  Program velocity_program;
  VAO velocity_vao;
  std::vector<glm::uvec3> vel_indices;
  std::vector<std::array<glm::vec3, 2>> vel_offsets;
};