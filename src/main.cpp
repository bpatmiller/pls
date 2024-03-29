#include <glad/glad.h>

#include "geom.h"
#include "gui.h"
#include <GLFW/glfw3.h>
#include <getopt.h>
#include <iostream>

void test_array3();

// ERROR CALLBACK
void glfw_error_callback(int error, const char *description) {
  fprintf(stderr, "Error %i: %s\n", error, description);
}

// KEY CALLBACKS
void KeyCallback(GLFWwindow *window, int key, int scancode, int action,
                 int mods) {

  GUI *gui = (GUI *)glfwGetWindowUserPointer(window);

  if (key == GLFW_KEY_Q) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
    return;
  } else if (key == GLFW_KEY_V && action == GLFW_RELEASE) {
    gui->dirty = true;
    gui->draw_velocity = !gui->draw_velocity;
  } else if (key == GLFW_KEY_B && action == GLFW_RELEASE) {
    gui->dirty = true;
    gui->draw_particles = !gui->draw_particles;
  } else if (key == GLFW_KEY_G && action == GLFW_RELEASE) {
    gui->dirty = true;
    gui->draw_grid = !gui->draw_grid;
  } else if (key == GLFW_KEY_F && action == GLFW_RELEASE) {
    gui->dirty = true;
    gui->display_phi = !gui->display_phi;
  } else if (key == GLFW_KEY_MINUS && action == GLFW_RELEASE) {
    // cycle through levelsets to draw
    gui->dirty = true;
    gui->levelset_to_draw =
        (gui->levelset_to_draw + 1) % gui->simulation.fluids.size();
  }
  if (action == GLFW_PRESS) {
    gui->keyHeld[key] = true;
  } else if (action == GLFW_RELEASE) {
    gui->keyHeld[key] = false;
  }
}

void MouseButtonCallback(GLFWwindow *window, int button, int action, int mods) {
  GUI *gui = (GUI *)glfwGetWindowUserPointer(window);
  if (button == GLFW_MOUSE_BUTTON_LEFT && action != GLFW_RELEASE) {
    gui->mouse_pressed = true;
  } else {
    gui->mouse_pressed = false;
  }
}

void MousePosCallback(GLFWwindow *window, double mouse_x, double mouse_y) {
  GUI *gui = (GUI *)glfwGetWindowUserPointer(window);
  gui->mouse_pos = glm::vec2(mouse_x, mouse_y);
}

int main(int argc, char *argv[]) {
  bool graphical = true;

  if (graphical) {
    // create window/init glfw
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      exit(EXIT_FAILURE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow *window = glfwCreateWindow(640, 480, "apic", NULL, NULL);
    if (!window) {
      throw std::runtime_error("glfwCreateWindow error");
    }
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSwapInterval(1);
    GUI gui(window);
    glfwSetWindowUserPointer(window, &gui);
    // key / mouse callbacks
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetCursorPosCallback(window, MousePosCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);

    gui.init(2.0f, 20, 20, 20);

    while (!glfwWindowShouldClose(window)) {
      gui.update();
      glfwSwapBuffers(window);
      glfwPollEvents();
    }
  } else {
    std::printf("Using headless mode\n");
    int n_dim = 50;
    Simulation simulation(2.0f, n_dim, n_dim, n_dim);
    simulation.add_two_walls();
    float time = 0.0f;
    float goal_time = 1.0f;
    float timestep = 0.1f;
    while (time < goal_time) {
      simulation.advance(timestep);
      time += timestep;
      simulation.export_mesh(time);
    }
  }
}

void test_array3() {
  std::cout << "TESTING ARRAY3 FUNCTIONS/OPERATORS\n";
  Array3f f(2, 2, 2);
  f.set(2);
  f.print();

  // assignment
  Array3f mis(2, 2, 2);
  mis.set(3);
  f = mis;
  // check that its not just a shallow copy
  mis.set(2);
  f.print();
  mis.set(3);

  // addition
  f.set(2);
  f = f + mis;
  // we should get an array of all 5
  f.print();

  f = f - 1.5;
  // we should get an array of all 3.5
  f.print();

  mis.set(2);
  f = f * mis;
  // array of all 7
  f.print();

  mis.set(0.5);
  f = f - mis;
  // array of 6.5
  f.print();

  mis.set(3);
  // 6.5 - (.5 * 2) = 5.5
  f = f - (mis - 1) * 0.5;
  f.print();

  // 11
  f = f * 2;
  f.print();
}