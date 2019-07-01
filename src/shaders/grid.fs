#version 430 core
out vec4 fragment_color;

in float p;

void main() {
  if (p < 0) {
    fragment_color = vec4(0.0, 1.0, 0.0, 1.0);
  } else {
    fragment_color = vec4(0.7, 0.1, p, 0.15);
  }
}