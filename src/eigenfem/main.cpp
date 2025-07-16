//
// main.cpp
//

#include <iostream>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"
#include "../../third-party/spectra-1.1.0/include/Spectra/GenEigsSolver.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseGenMatProd.h"

#include "materials.h"
#include "mesh.h"
#include "model.h"


int main()
{
    Material steel(7800., 2.1e11, 0.3);
    Material aluminium(2700., 7e10, 0.33);

    std::string mesh_path = "../data/bar.mesh";
    Material material = aluminium;
    Mesh mesh(mesh_path, material);
    mesh.import_gmsh_matlab();

    // std::cout << "Number of tetrahedra in mesh: " << mesh.n_elements << std::endl << std::endl;
    // std::cout << "Tetrahedra connectivity table:" << std::endl;
    // std::cout << mesh.table_tets << std::endl;

    std::vector<int> dirichlet_tags = {2};
    Model model(mesh, dirichlet_tags);
    model.create_dof_lists();
    std::cout << "Number of dirichlet dofs: " << model.dirichlet_dofs.size() << std::endl;
    //std::cout << "dirichlet dofs:" << std::endl;
    //for (size_t i = 0; i < model.dirichlet_dofs.size(); i++)
    //{
    //    std::cout << model.dirichlet_dofs[i] << std::endl;
    //}

    return 0;
}