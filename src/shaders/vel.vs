#version 430 core
layout(location = 0) in vec3 vertex_position;
layout(location = 1) in vec3 offset_position;
layout(location = 2) in vec3 offset_velocity;

uniform mat4 projection;
uniform mat4 view;

out vec3 v;

void main() {
  v = offset_velocity * 0.05;
  // manually handle all 4 vertices
  if (vertex_position.y == 0) {
    gl_Position =
        projection * view * vec4(offset_position + vertex_position, 1.0);
  } else {
    gl_Position = projection * view *
                  vec4(offset_position + offset_velocity * 0.01 +
                           (vertex_position - vec3(0, 1, 0)),
                       1.0);
  }
}
