//
// io.cpp
//

#include "io.h"


VTKwriter::VTKwriter(Mesh msh, Eigen::VectorXf vec_displacement)
{
    mesh = msh;
    U = vec_displacement;
}

void VTKwriter::add_U_to_mesh()
{
    deformed_mesh = mesh;
    Eigen::MatrixXf reshaped_U = U.reshaped(3, mesh.n_nodes).transpose();
}

void VTKwriter::write_deformed_mesh(std::string path_to_file, std::string filename_wo_extension, std::string trail)
{
    std::ofstream file;
    std::string full_path_to_file = path_to_file + filename_wo_extension + trail + ".vtk";
    file.open(full_path_to_file);

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << filename_wo_extension << ", Created by eigenfem" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    file.close();
}