#version 430 core
layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 position;
layout(location = 2) in float mass;
layout(location = 3) in vec3 velocity;
layout(location = 4) in int id;

uniform mat4 projection;
uniform mat4 view;

out vec3 v;
out flat int i;

void main() {
  v = velocity;
  i = id;
  gl_Position = projection * view * vec4(position + vertex_position, 1.0);
}
