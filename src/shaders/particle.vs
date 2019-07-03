#version 430 core
layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 position;
layout(location = 2) in float id;
layout(location = 3) in float phi;
layout(location = 4) in float radius;
layout(location = 5) in float valid;

uniform mat4 projection;
uniform mat4 view;

out float p;

void main() {
  p = phi;
  gl_Position = projection * view * vec4(position + vertex_position, 1.0);
}
