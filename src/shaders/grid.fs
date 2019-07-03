#version 430 core
out vec4 fragment_color;

uniform bool display_phi;
in float p;

void main() {
  // phi
  if (display_phi) {
    if (p < 0) {
      fragment_color = vec4(0.0, 1.0, 0.0, 1.0);
    } else if (p < 0.2) {
      fragment_color = vec4(0.7, 0.7, 0.2, 0.25);
    } else {
      fragment_color = vec4(0.7, 0.1, p, 0.15);
    }
  }
  // pressure
  else {
    if (p > 0) {
      fragment_color = vec4(0.6 + 5 * p, 0.3, 0.7, 0.5 + p);

    } else {
      fragment_color = vec4(0.1 + 2 * p, 0.3, 0.7, 0.5 + p);
    }
  }
}