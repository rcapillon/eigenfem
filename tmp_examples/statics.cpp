//
// statics.cpp
//

#include <iostream>
#include <ctime>

#include "solvers.h"
#include "io.h"


int main()
{
    time_t global_timer_start = time(nullptr); // Starts timer for whole code execution

    // Defines two basic metals as usable materials
    Material steel(7800., 2.1e11, 0.3); 
    Material aluminium(2700., 7e10, 0.33);

    std::string mesh_path = "../data/bar.mesh";
    Material material = aluminium;
    Mesh mesh(mesh_path, material); // Chosen material will be used for the whole domain
    mesh.import_gmsh_matlab(); // Constructs coordinates and connectivity tables

    std::vector<int> dirichlet_tags = {2}; // Specifies 2D physical group tags where 0-dirichlet conditions are imposed
    
    int surface_force_tag = 3; // Specifies 2D physical group tag where a surface force is applied

    // Defines the surface force vector
    Eigen::VectorXf vec_surface_force = Eigen::VectorXf::Zero(3);
    vec_surface_force(0) = 4e8;

    std::tuple<int, Eigen::VectorXf> tuple_surface_force = std::make_tuple(surface_force_tag, vec_surface_force);
    std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces;
    surf_forces.push_back(tuple_surface_force);

    int volume_force_tag = 1; // Specifies 3D physical group tag where a volume force is applied (can only be 1 for the moment)

    // Defines the volume force vector
    Eigen::VectorXf vec_volume_force = Eigen::VectorXf::Zero(3);
    // vec_volume_force(2) = 9.81 * 1e6;

    std::tuple<int, Eigen::VectorXf> tuple_volume_force = std::make_tuple(volume_force_tag, vec_volume_force);
    std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces;
    vol_forces.push_back(tuple_volume_force);

    Model model(mesh, dirichlet_tags, surf_forces, vol_forces);
    
    // LinearStaticsSolver can be used to solve a linear statics problem involving surface and volume forces
    LinearStaticsSolver solver(model);
    time_t solver_timer_start = time(nullptr); // Starts timer for solver execution
    solver.solve();
    time_t solver_timer_end = time(nullptr); // Ends timer for solver execution

    // Export deformed mesh to VTK format
    Eigen::VectorXf U = solver.U;
    VTKwriter vtk_writer(mesh, U);
    vtk_writer.add_U_to_mesh();

    // The first argument (path to folder) needs to end with "/"
    // The folders/subfolders also need to already be created.
    vtk_writer.write_deformed_mesh("../../eigenfem_vtk/", "statics");
    
    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << 
        std::endl << "Time spent solving the problem: " << solver_timer_end - solver_timer_start << " seconds." << std::endl;

    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}