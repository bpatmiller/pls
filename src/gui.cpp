#include "gui.h"

const glm::vec3 FORWARD(0, 0, -1);
const glm::vec3 SIDE(1, 0, 0);
const glm::vec3 UP(0, 1, 0);

void GUI::init(float lx_, int nx_, int ny_, int nz_) {
  // set up simulation
  simulation.init(lx_, nx_, ny_, nz_);

  // compile shaders
  grid_program = Program("src/shaders/grid.vs", "", "src/shaders/grid.fs", "");

  // scale the grids down to size
  float h = simulation.h;
  for (auto &v : box_vertices) {
    v *= h;
  }

  // add fluid-containing grid cells
  grid_offsets.resize(simulation.fluid_cell_count);
  int fcc = 0;
  for (int i = 1; i < simulation.liquid_phi.sx - 1; i++) {
    for (int j = 1; j < simulation.liquid_phi.sy - 1; j++) {
      for (int k = 1; k < simulation.liquid_phi.sz - 1; k++) {
        if (simulation.liquid_phi(i, j, k) < 0) {
          grid_offsets[fcc] = glm::vec4(h * i, h * j, h * k, 1.0f);
          fcc++;
        }
      }
    }
  }

  // set up grid VAO
  grid_vao.setLayout({3}, false);
  grid_vao.setLayout({4}, true);
  grid_vao.vb.set(box_vertices);
  grid_vao.ib.set(grid_offsets);

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

    // update grid vao
    grid_offsets.clear();
    grid_offsets.resize(simulation.fluid_cell_count);
    float h = simulation.h;
    int fcc = 0;
    for (int i = 1; i < simulation.liquid_phi.sx - 1; i++) {
      for (int j = 1; j < simulation.liquid_phi.sy - 1; j++) {
        for (int k = 1; k < simulation.liquid_phi.sz - 1; k++) {
          if (simulation.liquid_phi(i, j, k) < 0) {
            grid_offsets[fcc] = glm::vec4(h * i, h * j, h * k, 1.0f);
            fcc++;
          }
        }
      }
    }
    grid_vao.ib.set(grid_offsets, 0);
  }

  grid_program.use();
  grid_program.setMat4("projection", projection_matrix);
  grid_program.setMat4("view", view_matrix);
  grid_vao.bind();
  glDrawElementsInstanced(GL_TRIANGLES, box_indices.size() * 3, GL_UNSIGNED_INT,
                          box_indices.data(), grid_offsets.size());
}