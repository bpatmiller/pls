#version 430 core
layout(location = 0) in vec3 vertex_position;

uniform mat4 projection;
uniform mat4 view;

void main() { gl_Position = projection * view * vec4(vertex_position, 1.0); }
