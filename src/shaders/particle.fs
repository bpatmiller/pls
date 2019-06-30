#version 430 core
out vec4 fragment_color;

in vec3 v;
in flat int i;

void main() {
  // fragment_color =
  //     vec4(mix(vec3(0.1, 0.3, 0.7), vec3(0.05, 0.8, 0.5), length(v)), 1.0);
  fragment_color = vec4(sin(i), cos(i), i, 1.0);
}