#include "gui.h"

const glm::vec3 FORWARD(0, 0, -1);
const glm::vec3 SIDE(1, 0, 0);
const glm::vec3 UP(0, 1, 0);

void GUI::init(float lx_, int nx_, int ny_, int nz_) {
  // set up simulation
  simulation.init(lx_, nx_, ny_, nz_);
  // simulation.add_sphere_phi();
  simulation.add_pool_phi();

  // compile shaders
  grid_program = Program("src/shaders/grid.vs", "", "src/shaders/grid.fs", "");
  particle_program =
      Program("src/shaders/particle.vs", "", "src/shaders/particle.fs", "");
  velocity_program =
      Program("src/shaders/vel.vs", "", "src/shaders/vel.fs", "");

  // set up particle vao
  std::vector<glm::vec3> sphere_vertices;
  create_sphere(simulation.h * 0.1, sphere_vertices);
  particle_vao.setLayout({3}, false);
  particle_vao.setLayout({3, 1, 1, 1, 1}, true);
  particle_vao.vb.set(sphere_vertices);
  particle_vao.ib.set(simulation.fluids[0].particles);

  std::vector<glm::vec3> box_vertices = {
      {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f},
      {1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f},
      {1.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 1.0f}};
  // scale the grids down to size
  float h = simulation.h;
  for (auto &v : box_vertices) {
    v *= h;
  }

  // add grid cells
  grid_offsets.resize(simulation.fluids[0].phi.size);
  int fcc = 0;
  for (int i = 0; i < simulation.fluids[0].phi.sx; i++) {
    for (int j = 0; j < simulation.fluids[0].phi.sy; j++) {
      for (int k = 0; k < simulation.fluids[0].phi.sz; k++) {
        if (simulation.fluids[0].phi(i, j, k) < 0) {
          grid_offsets[fcc] = glm::vec4(h * i, h * j, h * k, 1.0f);
        } else {
          grid_offsets[fcc] = glm::vec4(h * i, h * j, h * k, 0.0f);
        }
        fcc++;
      }
    }
  }

  box_indices = {{0, 1, 2}, {1, 3, 2}, {4, 6, 5}, {5, 6, 7},
                 {0, 5, 1}, {0, 4, 5}, {2, 3, 7}, {2, 7, 6},
                 {3, 1, 5}, {3, 5, 7}, {0, 2, 6}, {0, 6, 4}};

  // set up grid VAO
  grid_vao.setLayout({3}, false);
  grid_vao.setLayout({4}, true);
  grid_vao.vb.set(box_vertices);
  grid_vao.ib.set(grid_offsets);

  // set up velocity field vao
  std::vector<glm::vec3> vel_vertices = {
      {-0.005f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.005f, 0.0f, 0.0f}};
  vel_indices = {{0, 1, 2}};
  vel_offsets.resize(simulation.fluids[0].phi.size);
  for (int i = 0; i < simulation.fluids[0].phi.sx; i++) {
    for (int j = 0; j < simulation.fluids[0].phi.sy; j++) {
      for (int k = 0; k < simulation.fluids[0].phi.sz; k++) {
        glm::vec3 p =
            glm::vec3(h * i + 0.5f * h, h * j + 0.5f * h, h * k + 0.5f * h);
        vel_offsets[i + (simulation.fluids[0].phi.sx * j) +
                    (simulation.fluids[0].phi.sx * simulation.fluids[0].phi.sy *
                     k)] = {p, glm::vec3(0)}; // simulation.trilerp_uvw(p)};
      }
    }
  }
  velocity_vao.setLayout({3}, false);
  velocity_vao.setLayout({3, 3}, true);
  velocity_vao.vb.set(vel_vertices);
  velocity_vao.ib.set(vel_offsets);

  // gl settings
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthFunc(GL_LESS);
}

void GUI::update_camera() {
  if (mouse_diff == glm::vec2(0, 0))
    return;
  mouse_diff *= -0.005f;

  pitch = std::max(-1.5f, std::min(1.5f, pitch + mouse_diff.y));
  yaw += mouse_diff.x;

  glm::quat qyaw = glm::angleAxis(yaw, UP);
  glm::quat qpitch = glm::angleAxis(pitch, SIDE);
  orientation = qyaw * qpitch;
}

void GUI::update(float t, bool force) {
  if (t < 0)
    t = timestep;
  // update camera vars
  glfwGetWindowSize(window, &window_dims.x, &window_dims.y);
  glViewport(0, 0, window_dims.x, window_dims.y);

  eye = focus + glm::vec3(glm::mat4_cast(orientation) *
                          glm::vec4(base_eye - focus, 1.0));
  view_matrix = glm::lookAt(eye, focus, UP);
  projection_matrix = glm::perspective(
      glm::radians(60.0f), ((float)window_dims.x) / window_dims.y, 0.01f, 30.f);

  // handle mouse movement
  bool first = (mouse_pos_prev == glm::vec2(-1, -1));
  if (mouse_pressed) {
    if (!first) {
      mouse_diff = mouse_pos - mouse_pos_prev;
      update_camera();
    }
  }
  mouse_pos_prev = mouse_pos;

  // handle keypress
  if (keyHeld[GLFW_KEY_W]) {
    if (base_eye.z - focus.z > 0.2f) {
      base_eye.z -= 0.1f;
    }
  }
  if (keyHeld[GLFW_KEY_S]) {
    base_eye.z += 0.1f;
  }

  // clear the renderer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // step the simulation, copy the new particle data
  if (keyHeld[GLFW_KEY_P] || force) {
    simulation.advance(t);
    dirty = true;
  }
  if (dirty) {
    // update grid vao
    float h = simulation.h;
    int fcc = 0;
    for (int i = 0; i < simulation.fluids[0].phi.sx; i++) {
      for (int j = 0; j < simulation.fluids[0].phi.sy; j++) {
        for (int k = 0; k < simulation.fluids[0].phi.sz; k++) {
          if (display_phi) {
            grid_offsets[fcc] =
                glm::vec4(h * i, h * j, h * k,
                          simulation.fluids[levelset_to_draw].phi(i, j, k));
          } else {
            grid_offsets[fcc] =
                glm::vec4(h * i, h * j, h * k, simulation.pressure(i, j, k));
          }
          fcc++;
        }
      }
    }
    grid_vao.ib.update(grid_offsets, 0);

    float offs = simulation.h * 0.5;
    // update vao
    for (int i = 0; i < simulation.fluids[0].phi.sx; i++) {
      for (int j = 0; j < simulation.fluids[0].phi.sy; j++) {
        for (int k = 0; k < simulation.fluids[0].phi.sz; k++) {
          glm::vec3 p =
              glm::vec3(simulation.h * i + offs, simulation.h * j + offs,
                        simulation.h * k + offs);

          vel_offsets[i + (simulation.fluids[0].phi.sx * j) +
                      (simulation.fluids[0].phi.sx *
                       simulation.fluids[0].phi.sy * k)][1] =
              simulation.trilerp_uvw(p);
        }
      }
    }
    velocity_vao.ib.update(vel_offsets, 0);

    particle_vao.setLayout({3, 1, 1, 1, 1}, true);
    particle_vao.ib.set(simulation.fluids[0].particles);
  }

  // draw velocity field
  if (draw_velocity) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    velocity_program.use();
    velocity_program.setMat4("projection", projection_matrix);
    velocity_program.setMat4("view", view_matrix);
    velocity_vao.bind();
    glDrawElementsInstanced(GL_TRIANGLES, vel_indices.size() * 3,
                            GL_UNSIGNED_INT, vel_indices.data(),
                            simulation.fluids[0].phi.size);
  }

  // draw grid
  if (draw_grid) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    grid_program.use();
    grid_program.setMat4("projection", projection_matrix);
    grid_program.setMat4("view", view_matrix);
    grid_program.setBool("display_phi", display_phi);
    grid_vao.bind();
    glDrawElementsInstanced(GL_TRIANGLES, box_indices.size() * 3,
                            GL_UNSIGNED_INT, box_indices.data(),
                            grid_offsets.size());
  }

  // draw particle
  if (draw_particles) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    particle_program.use();
    particle_program.setMat4("projection", projection_matrix);
    particle_program.setMat4("view", view_matrix);
    particle_vao.bind();
    glDrawElementsInstanced(GL_TRIANGLES, sphere_indices.size() * 3,
                            GL_UNSIGNED_INT, sphere_indices.data(),
                            simulation.fluids[0].particles.size());
  }

  dirty = false;
}

void GUI::create_sphere(float Radius, std::vector<glm::vec3> &s_vertices) {
  int Stacks = 4;
  int Slices = 4;
  s_vertices.clear();

  for (int i = 0; i <= Stacks; ++i) {
    float V = i / (float)Stacks;
    float phi = V * glm::pi<float>();
    for (int j = 0; j <= Slices; ++j) {
      float U = j / (float)Slices;
      float theta = U * (glm::pi<float>() * 2);
      float x = cosf(theta) * sinf(phi);
      float y = cosf(phi);
      float z = sinf(theta) * sinf(phi);
      s_vertices.push_back(glm::vec3(x, y, z) * Radius);
    }
  }

  for (int i = 0; i < Slices * Stacks + Slices; ++i) {
    sphere_indices.emplace_back(glm::uvec3(i, i + Slices + 1, i + Slices));
    sphere_indices.emplace_back(glm::uvec3(i + Slices + 1, i, i + 1));
  }
}