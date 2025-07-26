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

    // count and print of arguments fed to main
    // std::cout << "Number of arguments: " << argc << std::endl;
    // for (size_t i = 0; i < argc; i++)
    // {
    //     std::cout << "Argument " << i << ":" << std::endl;
    //     std::cout << argv[i] << std::endl << std::endl;
    // }

    std::string path_to_input_file = argv[1];

    // Parse input file
    InputParser ip(path_to_input_file);
    ip.parse_input_file();
    // for (size_t i = 0; i < ip.lines.size(); i++)
    // {
    //     std::cout << ip.lines[i] << std::endl;        
    // }

    std::cout << "mesh file path: " << ip.inputs.mesh_path << std::endl << std::endl;
    std::cout << "material rho: " << ip.inputs.material_rho << std::endl << std::endl;
    std::cout << "material young: " << ip.inputs.material_youngmodulus << std::endl << std::endl;
    std::cout << "material nu: " << ip.inputs.material_poissonratio << std::endl << std::endl;
    std::cout << "dirichlet tags: " << std::endl;
    for (size_t i = 0; i < ip.inputs.dirichlet_tags.size(); i++)
    {
        std::cout << ip.inputs.dirichlet_tags[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "volume force: " << ip.inputs.volume_force[0] << " , " << ip.inputs.volume_force[1] << " , " << ip.inputs.volume_force[2] << std::endl << std::endl;
    std::cout << "damping M: " << ip.inputs.damping_alpha_M << std::endl << std::endl;
    std::cout << "damping K: " << ip.inputs.damping_alpha_K << std::endl << std::endl;
    std::cout << "output name: " << ip.inputs.output_name << std::endl << std::endl;
    std::cout << "output path: " << ip.inputs.output_path << std::endl << std::endl;
    
    time_t global_timer_end = time(nullptr); // Ends timer for whole code execution

    // Prints timer values
    std::cout << "Global elapsed time: " << global_timer_end - global_timer_start << " seconds." << std::endl;

    return 0;
}