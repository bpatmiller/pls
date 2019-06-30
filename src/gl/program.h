#pragma once

#include <glad/glad.h>

#include <GLFW/glfw3.h>
#include <fstream>
#include <glm/glm.hpp>
#include <iostream>
#include <sstream>
#include <string>

class Program {
public:
  unsigned int ID;

  std::string readShader(const char *shaderPath) {
    std::string shaderCode;
    std::ifstream shaderFile;

    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      shaderFile.open(shaderPath);
      std::stringstream shaderStream;
      shaderStream << shaderFile.rdbuf();
      shaderFile.close();
      shaderCode = shaderStream.str();
    } catch (std::ifstream::failure &e) {
    }
    return shaderCode;
  }

  GLuint compile_shader(const char *shaderCode, GLint type,
                        std::string type_string, std::string fname) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &shaderCode, NULL);
    glCompileShader(s);
    checkCompileErrors(s, type_string, fname);
    return s;
  }

  Program() {}

  Program(const char *vertexPath, const char *geometryPath,
          const char *fragmentPath, const char *computePath) {
    std::string vertCode = readShader(vertexPath);
    std::string geomCode = readShader(geometryPath);
    std::string fragCode = readShader(fragmentPath);
    std::string compCode = readShader(computePath);

    GLuint vertex = 0, fragment = 0, geometry = 0, compute = 0;

    if (vertCode != "") {
      const char *vShaderCode = vertCode.c_str();
      vertex = compile_shader(vShaderCode, GL_VERTEX_SHADER, "VERTEX",
                              std::string(vertexPath));
    }
    if (geomCode != "") {
      const char *gShaderCode = geomCode.c_str();
      geometry = compile_shader(gShaderCode, GL_GEOMETRY_SHADER, "GEOMETRY",
                                std::string(geometryPath));
    }
    if (fragCode != "") {
      const char *fShaderCode = fragCode.c_str();
      fragment = compile_shader(fShaderCode, GL_FRAGMENT_SHADER, "FRAGMENT",
                                std::string(fragmentPath));
    }
    if (compCode != "") {
      const char *cShaderCode = compCode.c_str();
      compute = compile_shader(cShaderCode, GL_COMPUTE_SHADER, "COMPUTE",
                               std::string(computePath));
    }

    ID = glCreateProgram();

    if (vertCode != "") {
      glAttachShader(ID, vertex);
    }
    if (geomCode != "") {
      glAttachShader(ID, geometry);
    }
    if (fragCode != "") {
      glAttachShader(ID, fragment);
    }
    if (compCode != "") {
      glAttachShader(ID, compute);
    }

    glLinkProgram(ID);
    checkCompileErrors(ID, "PROGRAM");

    glDeleteShader(vertex);
    glDeleteShader(geometry);
    glDeleteShader(fragment);
    glDeleteShader(compute);
  }

  // activate the shader
  // ------------------------------------------------------------------------
  void use() { glUseProgram(ID); }
  // utility uniform functions
  // ------------------------------------------------------------------------
  void setBool(const std::string &name, bool value) const {
    glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value);
  }
  // ------------------------------------------------------------------------
  void setInt(const std::string &name, int value) const {
    glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
  }
  // ------------------------------------------------------------------------
  void setFloat(const std::string &name, float value) const {
    glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
  }
  // ------------------------------------------------------------------------
  void setMat4(const std::string &name, glm::mat4 value) const {
    glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE,
                       (const GLfloat *)&value[0][0]);
  }
  // ------------------------------------------------------------------------
  void setMat2(const std::string &name, glm::mat2 value) const {
    glUniformMatrix2fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE,
                       (const GLfloat *)&value[0][0]);
  }
  // ------------------------------------------------------------------------
  void setVec3(const std::string &name, glm::vec3 value) const {
    glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1,
                 (const GLfloat *)&value[0]);
  }
  // ------------------------------------------------------------------------
  void setVec2(const std::string &name, glm::vec2 value) const {
    glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1,
                 (const GLfloat *)&value[0]);
  }

private:
  // utility function for checking shader compilation/linking errors.
  // ------------------------------------------------------------------------
  void checkCompileErrors(unsigned int shader, std::string type,
                          std::string fname = "") {
    int success;
    char infoLog[1024];
    if (type != "PROGRAM") {
      glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
      if (!success) {
        glGetShaderInfoLog(shader, 1024, NULL, infoLog);
        std::cout
            << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n"
            << "FILENAME:: " << fname << std::endl
            << infoLog
            << "\n -- --------------------------------------------------- --"
            << std::endl;
      }
    } else {
      glGetProgramiv(shader, GL_LINK_STATUS, &success);
      if (!success) {
        glGetProgramInfoLog(shader, 1024, NULL, infoLog);
        std::cout
            << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n"
            << infoLog
            << "\n -- --------------------------------------------------- -- "
            << std::endl;
      }
    }
  }
};