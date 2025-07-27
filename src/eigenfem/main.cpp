//
// main.cpp
//

#include <iostream>
#include <ctime>

#include "solvers.h"
#include "io.h"


int main(int argc, char *argv[])
{
    time_t global_timer_start = time(nullptr); // Starts timer for whole code execution

    // Value used for pi
    const float PI = 3.141592653589793;

    std::string path_to_input_file = argv[1];

    // Parse input file
    InputParser ip(path_to_input_file);
    ip.parse_input_file();

    std::string mesh_path;
    if (ip.inputs.has_mesh)
    {
        mesh_path = ip.inputs.mesh_path;
    }
    else
    {
        mesh_path = "";
    }
    Material material;
    if (ip.inputs.has_material)
    {
        material = Material(ip.inputs.material_rho, 
            ip.inputs.material_youngmodulus, 
            ip.inputs.material_poissonratio);
    }
    else
    {
        material = Material();
    }

    Mesh mesh(mesh_path, material);

    std::vector<int> dirichlet_tags;
    if (ip.inputs.has_dirichlet)
    {
        dirichlet_tags = ip.inputs.dirichlet_tags;
    }
    else
    {
        dirichlet_tags = {};
    }
    std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces;
    std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces;
    if (ip.inputs.has_forces)
    {
        if (ip.inputs.has_volume_force)
        {
            Eigen::VectorXf vec_volume_force = Eigen::VectorXf::Zero(3);
            vec_volume_force(0) = ip.inputs.volume_force[0];
            vec_volume_force(1) = ip.inputs.volume_force[1];
            vec_volume_force(2) = ip.inputs.volume_force[2];
            std::tuple<int, Eigen::VectorXf> tuple_volume_force = std::make_tuple(1, vec_volume_force);
            vol_forces.push_back(tuple_volume_force);
        }
        else
        {
            vol_forces = {};
        }
        if (ip.inputs.has_surface_forces)
        {
            for (size_t i = 0; i < ip.inputs.tags_surface_forces.size(); i++)
            {
                Eigen::VectorXf vec_surface_force = Eigen::VectorXf::Zero(3);
                int surface_force_tag = ip.inputs.tags_surface_forces[i];
                vec_surface_force(0) = ip.inputs.surface_forces[i][0];
                vec_surface_force(1) = ip.inputs.surface_forces[i][1];
                vec_surface_force(2) = ip.inputs.surface_forces[i][2];
                std::tuple<int, Eigen::VectorXf> tuple_surface_force = std::make_tuple(surface_force_tag, vec_surface_force);
                surf_forces.push_back(tuple_surface_force);
            }
        }
        else
        {
            surf_forces = {};
        }
    }
    else
    {
        vol_forces = {};
        surf_forces = {};
    }
    float alpha_M;
    float alpha_K;
    if (ip.inputs.has_damping)
    {
        alpha_M = ip.inputs.damping_alpha_M;
        alpha_K = ip.inputs.damping_alpha_K;
    }
    else
    {
        alpha_M = 0.;
        alpha_K = 0.;
    }
    
    Model model(mesh, dirichlet_tags, surf_forces, vol_forces, alpha_M, alpha_K);

    if (ip.inputs.has_solver)
    {
        if (ip.inputs.solver_type.compare("MODAL") == 0)
        {
            int n_modes = ip.inputs.modal_n_modes;
            ModalSolver solver(model);
            solver.solve(n_modes);

            if (ip.inputs.has_output)
            {
                std::cout << "First " << n_modes << " smallest eigenfrequencies:" << std::endl;
                for (size_t i = 0; i < solver.vec_freqs.size(); i++)
                {
                    std::cout << solver.vec_freqs[i] << std::endl;
                }

                int plotted_mode_num = ip.inputs.plotted_mode_num;
                Eigen::VectorXf plotted_mode = solver.mat_modes(Eigen::all, plotted_mode_num);
                VTKwriter vtk_writer(mesh, plotted_mode);
                vtk_writer.add_U_to_mesh();
                vtk_writer.write_deformed_mesh(ip.inputs.output_path, ip.inputs.output_name);
            }
        }
        else if (ip.inputs.solver_type.compare("STATICS") == 0)
        {
            LinearStaticsSolver solver(model);
            solver.solve();

            if (ip.inputs.has_output)
            {
                Eigen::VectorXf U = solver.U;
                VTKwriter vtk_writer(mesh, U);
                vtk_writer.add_U_to_mesh();
                vtk_writer.write_deformed_mesh(ip.inputs.output_path, ip.inputs.output_name);
            }
            
        }
        else if (ip.inputs.solver_type.compare("FREQUENCYSWEEP") == 0)
        {
            FrequencySweepSolver solver(model);
            int n_freqs = ip.inputs.frequency_n_freq;
            float min_w = 2 * PI * ip.inputs.frequency_min_freq;
            float max_w = 2 * PI * ip.inputs.frequency_max_freq;
            float delta_w = (max_w - min_w) / (n_freqs - 1);
            float current_w = min_w - delta_w;
            std::vector<float> angular_freqs;
            for (size_t i = 0; i < n_freqs; i++)
            {
                current_w += delta_w;
                angular_freqs.push_back(current_w);
            }
            if (ip.inputs.solver_computes_rom)
            {
                int n_modes = ip.inputs.frequency_n_modes;
                solver.compute_rom(n_modes);
                solver.solve(angular_freqs);
            }
            else
            {
                DATio dat_io = DATio();
                Eigen::MatrixXf rom_basis = dat_io.load_dat(ip.inputs.frequency_path_to_basis) ;
                solver.load_rom(rom_basis);
                solver.solve(angular_freqs);
            }

            if (ip.inputs.has_output)
            {
                VTKwriter vtk_writer(mesh, solver.mat_U_modulus);
                vtk_writer.write_mesh_animation(ip.inputs.output_path, ip.inputs.output_name);
            }
        }
    }
    
    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}