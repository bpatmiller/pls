#pragma once

#include <glm/glm.hpp>

class Particle {
public:
  glm::vec3 position;
  glm::vec3 velocity;
  int id = 1;
  float phi;
  bool valid;

  Particle(glm::vec3 p) {
    position = p;
    velocity = glm::vec3(0, 0, 0);
  }

  Particle(glm::vec3 p, float ph) {
    position = p;
    phi = ph;
    valid = true;
  }
};