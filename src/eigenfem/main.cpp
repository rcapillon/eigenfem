//
// main.cpp
//

#include <iostream>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"
#include "../../third-party/spectra-1.1.0/include/Spectra/GenEigsSolver.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseGenMatProd.h"

#include "mesh.h"


int main()
{
    std::string mesh_path = "../data/bar.m";
    Mesh mesh(mesh_path);
    mesh.import_gmsh_matlab();

    std::cout << mesh.n_elements << std::endl << std::endl;
    std::cout << mesh.table_tets << std::endl;

    return 0;
}