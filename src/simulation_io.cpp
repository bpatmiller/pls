#include "marchingcubes.h"
#include "simulation.h"

void Simulation::export_mesh(float time) {
  Mesh m;
  for (int i = 0; i < (int)fluids.size(); i++) {
    m.vertices.clear();
    m.indices.clear();
    generate_mesh(fluids[i], m);
    m.save("fl" + std::to_string(i) + "_t" + std::to_string(time) + ".ply");
  }
}

void Simulation::generate_mesh(Fluid fluid, Mesh &mesh) {
  Array3f &phi = fluid.phi;

  std::vector<glm::vec3> positions(8);
  std::vector<float> values(8);
  int vert_count = 0;
  float offs = 0.5f * h;

  mesh.vertices.clear();
  mesh.indices.clear();

  int elements[8] = {0, 3, 4, 7, 1, 2, 5, 6};

  for (int i = 0; i < phi.sx - 1; i++) {
    for (int j = 0; j < phi.sy - 1; j++) {
      for (int k = 0; k < phi.sz - 1; k++) {
        positions.clear();
        values.clear();
        // pass in locations and values
        // of 8 neighboring sample points
        positions[elements[0]] =
            glm::vec3(i * h + offs, j * h + offs, k * h + offs);
        positions[elements[1]] =
            glm::vec3(i * h + offs, j * h + offs, (k + 1) * h + offs);
        positions[elements[2]] =
            glm::vec3(i * h + offs, (j + 1) * h + offs, k * h + offs);
        positions[elements[3]] =
            glm::vec3(i * h + offs, (j + 1) * h + offs, (k + 1) * h + offs);
        positions[elements[4]] =
            glm::vec3((i + 1) * h + offs, j * h + offs, k * h + offs);
        positions[elements[5]] =
            glm::vec3((i + 1) * h + offs, j * h + offs, (k + 1) * h + offs);
        positions[elements[6]] =
            glm::vec3((i + 1) * h + offs, (j + 1) * h + offs, k * h + offs);
        positions[elements[7]] = glm::vec3(
            (i + 1) * h + offs, (j + 1) * h + offs, (k + 1) * h + offs);

        values[elements[0]] = glm::clamp(phi(i, j, k), -0.5f, 0.5f);
        values[elements[1]] = glm::clamp(phi(i, j, k + 1), -0.5f, 0.5f);
        values[elements[2]] = glm::clamp(phi(i, j + 1, k), -0.5f, 0.5f);
        values[elements[3]] = glm::clamp(phi(i, j + 1, k + 1), -0.5f, 0.5f);
        values[elements[4]] = glm::clamp(phi(i + 1, j, k), -0.5f, 0.5f);
        values[elements[5]] = glm::clamp(phi(i + 1, j, k + 1), -0.5f, 0.5f);
        values[elements[6]] = glm::clamp(phi(i + 1, j + 1, k), -0.5f, 0.5f);
        values[elements[7]] = glm::clamp(phi(i + 1, j + 1, k + 1), -0.5f, 0.5f);

        polygonize(positions, values, vert_count, mesh.indices, mesh.vertices);
      }
    }
  }
}