#include <iostream>
#include <gmsh.h>
#include <vector>

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("Cilinder");

    double lc = 0.1;
    gmsh::model::geo::addPoint(0, 0, 0, lc, 1);  
    gmsh::model::geo::addPoint(1, 0, 0, lc, 2);
    gmsh::model::geo::addPoint(-1, 0, 0, lc, 3);

    gmsh::model::geo::addPoint(0, 0, 2, lc, 4);
    gmsh::model::geo::addPoint(1, 0, 2, lc, 5);
    gmsh::model::geo::addPoint(-1, 0, 2, lc, 6);

    gmsh::model::geo::addCircleArc(2, 1, 3, 1);
    gmsh::model::geo::addCircleArc(3, 1, 2, 2);

    gmsh::model::geo::addCircleArc(5, 4, 6, 3);
    gmsh::model::geo::addCircleArc(6, 4, 5, 4);

    gmsh::model::geo::addLine(2, 5, 5);
    gmsh::model::geo::addLine(3, 6, 6);

    gmsh::model::geo::addCurveLoop({1, 2}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);

    gmsh::model::geo::addCurveLoop({3, 4}, 2);
    gmsh::model::geo::addPlaneSurface({2}, 2);
    
    gmsh::model::geo::addCurveLoop({1, 6, -3, -5}, 3);
    gmsh::model::geo::addSurfaceFilling({3}, 3);
    gmsh::model::geo::addCurveLoop({2, 5, -4, -6}, 4);
    gmsh::model::geo::addSurfaceFilling({4}, 4);

    gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addVolume({1});

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("cilinder.msh");
    gmsh::finalize();
    return 0;
}