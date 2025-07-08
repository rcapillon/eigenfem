//
// main.cpp
//

#include <iostream>

#include "../../third-party/MshIO-0.0.1/include/mshio/mshio.h"
#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"
#include "../../third-party/spectra-1.1.0/include/Spectra/GenEigsSolver.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseGenMatProd.h"

#include "mesh.h"


int main()
{
    std::string mesh_path = "../data/bar.msh";
    Mesh mesh(mesh_path);
    mesh.import_gmsh();

    std::cout << mesh.nodes_table.table << std::endl;

    return 0;
}