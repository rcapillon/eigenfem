//
// example_frequency.cpp
//

#include <iostream>
#include <ctime>

#include "../src/eigenfem/solvers.h"
#include "../src/eigenfem/io.h"


int main()
{
    time_t global_timer_start = time(nullptr); // Starts timer for whole code execution

    // Value used for pi
    const float PI = 3.141592653589793;

    // Defines two basic metals as usable materials
    Material steel(7800., 2.1e11, 0.3); 
    Material aluminium(2700., 7e10, 0.33);

    std::string mesh_path = "../examples/data/plate.mesh";
    Material material = aluminium;
    Mesh mesh(mesh_path, material); // Chosen material will be used for the whole domain
    mesh.import_gmsh_matlab(); // Constructs coordinates and connectivity tables

    std::vector<int> dirichlet_tags = {2}; // Specifies 2D physical group tags where 0-dirichlet conditions are imposed
    
    int surface_force_tag = 3; // Specifies 2D physical group tag where a surface force is applied

    // Defines the surface force vector
    Eigen::VectorXf vec_surface_force = Eigen::VectorXf::Zero(3);
    vec_surface_force(0) = 1.25e10;

    std::tuple<int, Eigen::VectorXf> tuple_surface_force = std::make_tuple(surface_force_tag, vec_surface_force);
    std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces;
    surf_forces.push_back(tuple_surface_force);

    // Rayleigh damping coefficients
    float alpha_M = 1e-6;
    float alpha_K = 1e-6;

    Model model(mesh, dirichlet_tags, surf_forces, {}, alpha_M, alpha_K);
    
    // FrequencySweepSolver can be used to solve a frequency-domain dynamic problem 
    // using a reduced-order model and involving surface and volume forces
    FrequencySweepSolver solver(model);
    int n_modes = 10;
    int n_freqs = 300;
    float min_w = 2 * PI * 1200;
    float max_w = 2 * PI * 2200;
    float delta_w = (max_w - min_w) / (n_freqs - 1);
    float current_w = min_w - delta_w;
    std::vector<float> angular_freqs;
    for (size_t i = 0; i < n_freqs; i++)
    {
        current_w += delta_w;
        angular_freqs.push_back(current_w);
    }
    
    time_t rom_timer_start = time(nullptr); // Starts timer for reduced-order model construction
    std::cout << "Computing ROM..." << std::endl;
    solver.compute_rom(n_modes);
    time_t rom_timer_end = time(nullptr); // Ends timer for reduced-order model construction
    time_t solver_timer_start = time(nullptr); // Starts timer for solver execution
    std::cout << "Solving for all frequencies..." << std::endl;
    solver.solve(angular_freqs);
    time_t solver_timer_end = time(nullptr); // Ends timer for solver execution

    // Export deformed meshes of frequency sweep to VTK format
    std::cout << "Exporting frequency sweep to group of VTK files..." << std::endl;
    VTKwriter vtk_writer(mesh, solver.mat_U_modulus);

    // The first argument (path to folder) needs to end with "/"
    // The folders/subfolders also need to already be created.
    vtk_writer.write_mesh_animation("../examples/results/example_frequency_coarse/", "frequencysweep");
    
    std::cout << "Done." << std::endl;

    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << 
    std::endl << "Time spent constructing the reduced-order model: " << rom_timer_end - rom_timer_start << " seconds." << std::endl;
    std::cout << "Time spent solving the problem: " << solver_timer_end - solver_timer_start << " seconds." << std::endl;
    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}