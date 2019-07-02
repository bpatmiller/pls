# project defn
cmake_minimum_required(VERSION 3.10)
set(PROJECT_NAME pls)
project(${PROJECT_NAME})

# cmake options
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# cpp version / settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-fopenmp -Wall -Wextra -Werror -fmax-errors=1 -Wno-unused-parameter -Wno-deprecated-copy")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# source / lib directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
aux_source_directory(${SRC_DIR} SOURCES)
add_executable(${PROJECT_NAME} ${SOURCES})

# GLFW
find_package(PkgConfig REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
include_directories(${GLFW_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${GLFW_STATIC_LIBRARIES})

# GLM
FIND_PACKAGE(glm 0.9.8.0 QUIET)
INCLUDE_DIRECTORIES(${glm_INCLUDE_DIR})
add_definitions(-DGLM_ENABLE_EXPERIMENTAL -DGLM_FORCE_SIZE_FUNC=1 -DGLM_FORCE_RADIANS=1)

# GLAD
set(glad_dir ${CMAKE_SOURCE_DIR}/lib/glad)
add_library(glad ${glad_dir}/src/glad.c)
target_include_directories(glad PRIVATE ${glad_dir}/include)
target_include_directories(${PROJECT_NAME} PRIVATE ${glad_dir}/include)
target_link_libraries(${PROJECT_NAME} glad dl)