//
// io.cpp
//

#include "io.h"


VTKwriter::VTKwriter(Mesh msh, Eigen::VectorXf vec_displacement)
{
    mesh = msh;
    U = vec_displacement;
}

void VTKwriter::write_deformed_mesh(std::string path_to_file)
{
    // std::ofstream myfile;
    // myfile.open("example.txt");
    // myfile << "Writing this to a file.\n";
    // myfile.close();
}