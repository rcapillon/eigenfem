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
#include "utils.h"


int main()
{
    Material steel(7800., 2.1e11, 0.3);
    Material aluminium(2700., 7e10, 0.33);

    std::string mesh_path = "../data/bar.mesh";
    Material material = aluminium;
    Mesh mesh(mesh_path, material);
    mesh.import_gmsh_matlab();

    std::vector<int> dirichlet_tags = {2, 3};
    Model model(mesh, dirichlet_tags);
    model.create_dof_lists();
    std::cout << "Number of dirichlet dofs: " << model.dirichlet_dofs.size() << std::endl;
    std::cout << "Number of free dofs: " << model.free_dofs.size() << std::endl;

    const int n = 1000;
    Eigen::SparseMatrix<float> mat1(n, n);
    mat1.reserve(Eigen::VectorXi::Constant(n, 3));
    for(int i = 0; i < n; i++)
    {
        mat1.insert(i, i) = 1.0;
        if(i > 0)
            mat1.insert(i - 1, i) = 3.0;
        if(i < n - 1)
            mat1.insert(i + 1, i) = 2.0;
    }

    std::vector<int> inds = {0, 1, 2, 3, 4};
    SpMat mat2 = double_slice_spmat(mat1, inds, inds);

    std::cout << "sparse slice test : " << std::endl;
    std::cout << Eigen::MatrixXf(mat2) << std::endl;

    return 0;
}