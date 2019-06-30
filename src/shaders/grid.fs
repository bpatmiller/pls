#version 430 core
out vec4 fragment_color;

in float p;

void main() {
  // signed distance
  // if (p <= 0) {
  //   fragment_color = vec4(1.0);
  // } else {
  //   fragment_color = vec4(0.0, 0.5, (3 - p) * 0.25, clamp(3 - p, 0.15, 0.5));
  // }
  // pressure
  // if (p > 0) {
  //   fragment_color = vec4(p * .01, 0.1, 1 - (p * 0.01), 0.8);
  // } else {
  //   fragment_color = vec4(0.1);
  //}
  // divergence
  // if (p == 0) {
  //   //    fragment_color = vec4(0.8,0.8,0.8,0.05);
  //   fragment_color = vec4(0);
  // } else {
  //   fragment_color = vec4(sin(p), 0.2, 0.2, 0.9);
  // }

  // // type
  if (p == 0) {
    fragment_color = vec4(0.1, 0.8, 0.1, 0.4);
  } else if (p == 1) {
    fragment_color = vec4(0.2, 0.6, 0.6, 0.6);
  } else {
    fragment_color = vec4(0.7, 0.1, 0.4, 0.5);
  }
}