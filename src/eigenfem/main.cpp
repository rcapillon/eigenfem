//
// main.cpp
//

#include <iostream>

#include "solvers.h"


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
    
    ModalSolver solver(model);
    int n_modes = 5;
    solver.solve(n_modes);

    std::cout << "First " << n_modes << " smallest eigenfrequencies:" << std::endl;
    for (size_t i = 0; i < solver.vec_freqs.size(); i++)
    {
        std::cout << solver.vec_freqs[i] << std::endl;
    }
    
    return 0;
}