#version 430 core
out vec4 fragment_color;

uniform bool display_phi;
in float p;

void main() {
  // phi
  if (display_phi) {
    if (p < 0) {
      fragment_color = vec4(0.1, 1.0, 0 - 5 * p, 0.5);
    } else if (p < 0.2) {
      fragment_color = vec4(0.7, 0.7 - 2 * p, 0.2, 0.25);
    } else {
      fragment_color = vec4(0.7, 0.1, p, 0.15);
    }
  }
  // pressure
  else {
    if (abs(p) < 0.00001) {
      fragment_color = vec4(0.1 + 2 * p, 0.3, 0.7, 0.25);
    } else if (p > 0) {
      fragment_color = vec4(0.6 + p, 0.3, 0.7, 0.5 + p);

    } else {
      fragment_color = vec4(0.1 + p, 0.7, 0.7, 0.5 + abs(p));
    }
  }
}