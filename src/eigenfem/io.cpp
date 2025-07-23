//
// io.cpp
//

#include "io.h"


std::vector<float> read_matrix_row_from_line(std::string line)
{
	std::vector<float> row;
	size_t prev = 0;
	size_t pos;
	std::string substr;
    while ((pos = line.find_first_of(" ,", prev)) != std::string::npos)
    {
        if (pos > prev)
            row.push_back(std::stof(line.substr(prev, pos-prev)));
        prev = pos + 1;
    }
    if (prev < line.length())
        row.push_back(std::stof(line.substr(prev, std::string::npos)));

	return row;
}

VTKwriter::VTKwriter(Mesh msh, Eigen::VectorXf vec_displacement)
{
    mesh = msh;
    U = vec_displacement;
    deformed_mesh = mesh;
}

VTKwriter::VTKwriter(Mesh msh, Eigen::MatrixXf mat_displacement)
{
    mesh = msh;
    mat_U = mat_displacement;
    deformed_mesh = mesh;
}

void VTKwriter::add_U_to_mesh()
{
    Eigen::MatrixXf reshaped_U = U.reshaped(3, mesh.n_nodes).transpose();
    deformed_mesh.table_nodes += reshaped_U;
}

void VTKwriter::add_U_to_mesh(int idx_U)
{
    Eigen::MatrixXf reshaped_U = mat_U(Eigen::all, idx_U).reshaped(3, mesh.n_nodes).transpose();
    deformed_mesh.table_nodes += reshaped_U;
}

void VTKwriter::reset_deformed_mesh()
{
    deformed_mesh = mesh;
}

void VTKwriter::write_deformed_mesh(std::string path_to_folder, std::string filename_wo_extension, std::string trail)
{
    std::ofstream file;
    std::string full_path_to_file = path_to_folder + filename_wo_extension + trail + ".vtk";
    file.open(full_path_to_file);

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << filename_wo_extension << ", Created by eigenfem" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << deformed_mesh.n_nodes << " float" << std::endl;

    // POINTS
    for (size_t i = 0; i < deformed_mesh.table_nodes.rows(); i++)
    {
        for (size_t j = 0; j < deformed_mesh.table_nodes.cols(); j++)
        {
            file << deformed_mesh.table_nodes(i, j) << " ";
        }
        file << std::endl;
    }
    file << std::endl;

    // CELLS
    int n_tets = deformed_mesh.table_tets.rows();
    int n_tris = 0;
    for (size_t i = 0; i < deformed_mesh.tables_tris.size(); i++)
    {
        n_tris += deformed_mesh.tables_tris[i].rows();
    }
    int n_cells = n_tets + n_tris;
    int next_number = n_tets * 5 + n_tris * 4;

    file << "CELLS " << n_cells << " " << next_number << std::endl;
        // TRI CELLS
    for (size_t i = 0; i < deformed_mesh.tables_tris.size(); i++)
    {
        for (size_t j = 0; j < deformed_mesh.tables_tris[i].rows(); j++)
        {
            file << "3 ";
            for (size_t k = 0; k < deformed_mesh.tables_tris[i].cols(); k++)
            {
                file << deformed_mesh.tables_tris[i](j, k) << " ";
            }
            file << std::endl;
        }   
    }
        // TET CELLS
    for (size_t i = 0; i < deformed_mesh.table_tets.rows(); i++)
    {
        file << "4 ";
        for (size_t j = 0; j < deformed_mesh.table_tets.cols(); j++)
        {
            file << deformed_mesh.table_tets(i, j) << " ";
        }
        file << std::endl;
    }
    file << std::endl;
    
    // CELL_TYPES
    file << "CELL_TYPES " << n_cells << std::endl;
        // TRI CELLS
    for (size_t i = 0; i < deformed_mesh.tables_tris.size(); i++)
    {
        for (size_t j = 0; j < deformed_mesh.tables_tris[i].rows(); j++)
        {
            file << "5" << std::endl;
        }   
    }
        // TET CELLS
    for (size_t i = 0; i < deformed_mesh.table_tets.rows(); i++)
    {
        file << "10" << std::endl;
    }
    file << std::endl;
    
    // CELL_DATA
    file << "CELL_DATA " << n_cells << std::endl;
    file << "SCALARS CellEntityIds int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < deformed_mesh.tables_tris.size(); i++)
    {
        for (size_t j = 0; j < deformed_mesh.tables_tris[i].rows(); j++)
        {
            file << deformed_mesh.tris_tags[i] << std::endl;
        }   
    }
        // TET CELLS
    for (size_t i = 0; i < deformed_mesh.table_tets.rows(); i++)
    {
        file << "1" << std::endl;
    }
    file << std::endl;

    file.close();
}

void VTKwriter::write_mesh_animation(std::string path_to_folder, std::string filename_wo_extension)
{
    int n_freqs = mat_U.cols();
    std::string str_n_freqs = std::to_string(n_freqs);
    size_t max_digits = str_n_freqs.length();
    for (size_t i = 0; i < n_freqs; i++)
    {
        std::string trail = std::to_string(i + 1);
        size_t n_leading_zeros = max_digits - std::min(max_digits, trail.length());
        std::string complete_trail = "_" + std::string(n_leading_zeros, '0').append(trail);

        reset_deformed_mesh();
        add_U_to_mesh(i);
        write_deformed_mesh(path_to_folder, filename_wo_extension, complete_trail);
    }   
}

DATio::DATio(Eigen::MatrixXf mat)
{
    matrix = mat;
}

void DATio::write_matrix(std::string path_to_file)
{
    std::ofstream file;
    file.open(path_to_file);

    file << "ROWS" << std::endl;
    file << matrix.rows() << std::endl;
    file << "COLS" << std::endl;
    file << matrix.cols() << std::endl << std::endl;

    for (size_t i = 0; i < matrix.rows(); i++)
    {
        for (size_t j = 0; j < matrix.cols(); j++)
        {
            file << matrix(i, j) << " , ";
        }
        file << std::endl;
    }
    file << "END" << std::endl;
}

Eigen::MatrixXf DATio::load_dat(std::string path_to_file)
{
    std::ifstream file(path_to_file);
    std::string current_line;

    std::getline(file, current_line);
    std::getline(file, current_line);
    int n_rows = std::stoi(current_line);
    std::getline(file, current_line);
    std::getline(file, current_line);
    int n_cols = std::stoi(current_line);
    std::getline(file, current_line);

    Eigen::MatrixXf mat(n_rows, n_cols);

    for (size_t i = 0; i < n_rows; i++)
    {
        std::getline(file, current_line);
        std::vector<float> row = read_matrix_row_from_line(current_line);
        for (size_t j = 0; j < row.size(); j++)
        {
            mat(i, j) = row[j];
        }
    }
    
    return mat;
}

InputParser::InputParser(std::string path_to_input_file)
{

}