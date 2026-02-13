#include <iostream>
#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("stl");

    gmsh::merge("../livesey.stl"); 

    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;

    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);

    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    std::vector<int> surfaceTags;
    for (auto const& e : surfaces) surfaceTags.push_back(e.second);
 
    int l_Tag = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({l_Tag});

    gmsh::model::geo::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", 2.0);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 2.0);

    gmsh::model::mesh::generate(3);
    gmsh::write("stl.msh");
    gmsh::finalize();
    return 0;
}