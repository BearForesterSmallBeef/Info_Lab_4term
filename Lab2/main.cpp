#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

class CalcNode {
    friend class CalcMesh;
protected:
    double x, y, z;
    double smth;
    double vx, vy, vz;
    double x0, y0, z0;
public:
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0), x0(0.0), y0(0.0), z0(0.0)
    {
    }
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz)
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz), x0(x), y0(y), z0(z) {}

    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

class Element {
    friend class CalcMesh;
protected:
    unsigned long nodesIds[4]; 
};

class CalcMesh {
protected:
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {
        nodes.reserve(nodesCoords.size() / 3);
        for (size_t i = 0; i < nodesCoords.size(); i += 3) {
            double x = nodesCoords[i];
            double y = nodesCoords[i + 1];
            double z = nodesCoords[i + 2];
            
            // начальное скалярное поле
            double smth = sqrt(x*x + y*y + z*z); 
            
            nodes.emplace_back(x, y, z, smth, 0.0, 0.0, 0.0);
        }

        elements.resize(tetrsPoints.size() / 4);
        for (size_t i = 0; i < elements.size(); i++) {
            elements[i].nodesIds[0] = tetrsPoints[i * 4 + 0] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
        }
    }

    void doTimeStep(double tau, double time) {
        double jumpHeight = 50;
        double jumpSpeed = 25.0;
            
        double verticalPos = abs(sin(time * jumpSpeed)) * jumpHeight;
            
        double squish = 1.0 + 0.2 * sin(2*time * jumpSpeed - M_PI/2);

        for (auto& node : nodes) {
            node.x = node.x0 * (1.0 / sqrt(squish));
            node.y = node.y0 * (1.0 / sqrt(squish));
            node.z = node.z0 * squish;

            node.z += verticalPos;

            node.vz = -1.0 * jumpSpeed * jumpHeight * cos(2*time * jumpSpeed);
            node.vx = 0;
            node.vy = 0;

            // node.smth = sin(node.z * 0.01 - time * 10.0);
            node.smth = sin(node.z * 0.1 - time * 4.0);
        }
    }
// Я ЗАКОНЧИЛ ЗДЕСЬ!!!!!! 13 февраля 
    void snapshot(unsigned int snap_number) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();
        
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("my_scalar_field");
        
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for (const auto& node : nodes) {
            dumpPoints->InsertNextPoint(node.x, node.y, node.z);
            double _vel[3] = {node.vx, node.vy, node.vz};
            vel->InsertNextTuple(_vel);
            smth->InsertNextValue(node.smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);
       // чото байду я написал, мб там не спроста unsigned 
        for (const auto& elem : elements) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            for (int i = 0; i < 4; i++) 
                tetra->GetPointIds()->SetId(i, elem.nodesIds[i]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "lab_step_" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("main_lab");

    gmsh::merge("../Om_Nom.stl"); 

    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();
    
    //gmsh::option::setNumber("Mesh.MeshSizeMin", 2); 
    //gmsh::option::setNumber("Mesh.MeshSizeMax", 2);

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "2");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    // я не понял логику лабы, зачем мы так бережно относимся к четверке - надо спросить!!!
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] == 4) {
            tetrsNodesTags = &elementNodeTags[i];
            break;
        }
    }

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);
    
    gmsh::finalize();

    double tau = 0.01;
    double total_time = 0;
    
    cout << "Starting simulation" << endl;
    for(int step = 0; step < 100; step++) {
        mesh.doTimeStep(tau, total_time);
        mesh.snapshot(step);
        total_time += tau;
        if (step % 10 == 0) cout << "Step " << step << " done" << endl;
    }
    cout << "Чики-брики." << endl;

    return 0;
}