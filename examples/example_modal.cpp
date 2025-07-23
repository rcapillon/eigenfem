//
// example_modal.cpp
//

#include <iostream>
#include <ctime>

#include "../src/eigenfem/solvers.h"
#include "../src/eigenfem/io.h"


int main()
{
    time_t global_timer_start = time(nullptr); // Starts timer for whole code execution

    // Defines two basic metals as usable materials
    Material steel(7800., 2.1e11, 0.3); 
    Material aluminium(2700., 7e10, 0.33);

    std::string mesh_path = "../examples/data/bar.mesh";
    Material material = aluminium;
    Mesh mesh(mesh_path, material); // Chosen material will be used for the whole domain
    mesh.import_gmsh_matlab(); // Constructs coordinates and connectivity tables

    std::vector<int> dirichlet_tags = {2, 3}; // Specifies 2D physical group tags where 0-dirichlet conditions are imposed
    Model model(mesh, dirichlet_tags);
    
    ModalSolver solver(model); // ModalSolver can be used to compute the first N smallest eigenfrequencies and eigenvectors
    int n_modes = 5;
    time_t solver_timer_start = time(nullptr); // Starts timer for solver execution
    std::cout << "Calculating elastic modes..." << std::endl;
    solver.solve(n_modes);
    time_t solver_timer_end = time(nullptr); // Ends timer for solver execution

    // Print eigenfrequencies calculated by the solver
    std::cout << "First " << n_modes << " smallest eigenfrequencies:" << std::endl;
    for (size_t i = 0; i < solver.vec_freqs.size(); i++)
    {
        std::cout << solver.vec_freqs[i] << std::endl;
    }

    // Export a single mode shape to VTK format
    int plotted_mode_num = 0;
    std::cout << "Exporting mode " << plotted_mode_num << " to VTK..." << std::endl;
    Eigen::VectorXf plotted_mode = solver.mat_modes(Eigen::all, plotted_mode_num);
    VTKwriter vtk_writer(mesh, plotted_mode);
    vtk_writer.add_U_to_mesh();

    // The first argument (path to folder) needs to end with "/"
    // The folders/subfolders also need to already be created.
    vtk_writer.write_deformed_mesh("../examples/results/example_modal/", "mode0_shape");
    
    std::cout << "Done." << std::endl;

    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << 
    std::endl << "Time spent solving the problem: " << solver_timer_end - solver_timer_start << " seconds." << std::endl;
    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}