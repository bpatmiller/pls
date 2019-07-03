#pragma once

#include <glm/glm.hpp>

class Particle {
public:
  glm::vec3 position;
  int id;
  float starting_phi;
  float radius;
  bool valid;

  Particle(glm::vec3 p, float phi, float rad) {
    position = p;
    starting_phi = phi;
    radius = rad;
    valid = true;
  }
};