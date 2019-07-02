#include "gui.h"

const glm::vec3 FORWARD(0, 0, -1);
const glm::vec3 SIDE(1, 0, 0);
const glm::vec3 UP(0, 1, 0);

void GUI::init(float lx_, int nx_, int ny_, int nz_) {
  // set up simulation
  simulation.init(lx_, nx_, ny_, nz_);

  // compile shaders
  grid_program = Program("src/shaders/grid.vs", "", "src/shaders/grid.fs", "");

  std::vector<glm::vec3> box_vertices = {
      {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f},
      {1.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 1.0f, 1.0f},
      {1.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 1.0f}};
  // scale the grids down to size
  float h = simulation.h;
  for (auto &v : box_vertices) {
    v *= h;
  }

  // add fluid-containing grid cells
  grid_offsets.resize(simulation.liquid_phi.size);
  int fcc = 0;
  for (int i = 0; i < simulation.liquid_phi.sx; i++) {
    for (int j = 0; j < simulation.liquid_phi.sy; j++) {
      for (int k = 0; k < simulation.liquid_phi.sz; k++) {
        if (simulation.liquid_phi(i, j, k) < 0) {
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
    float h = simulation.h;
    int fcc = 0;
    for (int i = 0; i < simulation.liquid_phi.sx; i++) {
      for (int j = 0; j < simulation.liquid_phi.sy; j++) {
        for (int k = 0; k < simulation.liquid_phi.sz; k++) {
          grid_offsets[fcc] =
              glm::vec4(h * i, h * j, h * k, simulation.liquid_phi(i, j, k));
          fcc++;
        }
      }
    }
    grid_vao.ib.update(grid_offsets, 0);
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  grid_program.use();
  grid_program.setMat4("projection", projection_matrix);
  grid_program.setMat4("view", view_matrix);
  grid_vao.bind();
  glDrawElementsInstanced(GL_TRIANGLES, box_indices.size() * 3, GL_UNSIGNED_INT,
                          box_indices.data(), grid_offsets.size());
}