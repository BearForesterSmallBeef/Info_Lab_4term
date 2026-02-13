#include <iostream>
#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("Torus");

    int v1 = gmsh::model::occ::addTorus(0, 0, 0, 0.4, 0.2);
    
    int v2 = gmsh::model::occ::addTorus(0, 0, 0, 0.4, 0.15);

    std::vector<std::pair<int, int>> out;
    std::vector<std::vector<std::pair<int, int>>> outMap;
    gmsh::model::occ::cut({{3, v1}}, {{3, v2}}, out, outMap);

    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", 0.05);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.05);

    gmsh::model::mesh::generate(3);
    gmsh::write("torus.msh");
    gmsh::finalize();
    return 0;
}
