#pragma once

#include <glad/glad.h>

#include <GLFW/glfw3.h>
#include <glm/gtx/io.hpp>
#include <vector>

struct VBO {
  GLuint id = 0;

  VBO() { glGenBuffers(1, &id); }

  ~VBO() { glDeleteBuffers(1, &id); }

  void bind() { glBindBuffer(GL_ARRAY_BUFFER, id); }

  static void unbind() { glBindBuffer(GL_ARRAY_BUFFER, 0); }

  template <typename T>
  VBO &set(const std::vector<T> &vertices, GLenum usage = GL_STATIC_DRAW) {
    bind();
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(T), vertices.data(),
                 usage);
    unbind();
    return *this;
  }

  template <typename T>
  VBO &update(const std::vector<T> &vertices, size_t index) {
    bind();
    glBufferSubData(GL_ARRAY_BUFFER, index, vertices.size() * sizeof(T),
                    vertices.data());
    unbind();
    return *this;
  }
};

struct VAO {
  unsigned int vao;
  VBO vb;
  VBO ib;
  int attCt = 0;

  VAO() { glGenVertexArrays(1, &vao); }

  ~VAO() { glDeleteVertexArrays(1, &vao); }

  void bind() { glBindVertexArray(vao); }

  static void unbind() { glBindVertexArray(0); }

  void setLayout(std::initializer_list<uint> sizes, bool instanced) {

    bind();
    if (instanced) {
      ib.bind();
    } else {
      vb.bind();
    }

    int stride = 0;
    for (size_t size : sizes) {
      stride += size;
    }

    size_t offset = 0;
    for (int size : sizes) {
      glEnableVertexAttribArray(attCt);
      glVertexAttribPointer(attCt, size, GL_FLOAT, GL_FALSE,
                            stride * sizeof(float), (void *)offset);
      if (instanced) {
        glVertexAttribDivisor(attCt, 1);
      }
      offset += size * sizeof(float);
      attCt++;
    }

    VBO::unbind();
    unbind();
  }
};