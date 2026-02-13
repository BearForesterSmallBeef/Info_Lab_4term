#include <iostream>
#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("Circule");

    double lc = 0.1;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1); 
    gmsh::model::geo::addPoint(1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(-1, 0, 0, lc, 3);

    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 2, 2);

    gmsh::model::geo::addCurveLoop({1, 2}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);
    gmsh::write("circule.msh");
    gmsh::finalize();
    return 0;
}