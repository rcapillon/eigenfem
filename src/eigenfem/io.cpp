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

std::vector<int> read_ints_from_line(std::string line)
{
	std::vector<int> row;
	size_t prev = 0;
	size_t pos;
	std::string substr;
    while ((pos = line.find_first_of(" ,", prev)) != std::string::npos)
    {
        if (pos > prev)
            row.push_back(std::stoi(line.substr(prev, pos-prev)));
        prev = pos + 1;
    }
    if (prev < line.length())
        row.push_back(std::stoi(line.substr(prev, std::string::npos)));

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
    path_to_file = path_to_input_file;
}

void InputParser::get_lines(std::ifstream &file, std::string &line)
{
    while (line.compare("# END FILE") != 0)
    {
        std::getline(file, line);
        lines.push_back(line);
    }
    
}

void InputParser::parse_mesh()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# MESH") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## PATH") != 0 && current_line.compare("# END MESH") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## PATH") == 0)
    {
        inputs.has_mesh = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.mesh_path = current_line; 
    }
}

void InputParser::parse_material()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# MATERIAL") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## MASS DENSITY") != 0 && current_line.compare("# END MATERIAL") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## MASS DENSITY") == 0)
    {
        inputs.has_material = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.material_rho = std::stof(current_line);
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("## YOUNG MODULUS") != 0 && current_line.compare("# END MATERIAL") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## YOUNG MODULUS") == 0)
    {
        inputs.has_material = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.material_youngmodulus = std::stof(current_line); 
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("## POISSON RATIO") != 0 && current_line.compare("# END MATERIAL") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## POISSON RATIO") == 0)
    {
        inputs.has_material = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.material_poissonratio = std::stof(current_line); 
    }
}

void InputParser::parse_dirichlet()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# DIRICHLET") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## TAGS") != 0 && current_line.compare("# END DIRICHLET") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## TAGS") == 0)
    {
        inputs.has_dirichlet = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.dirichlet_tags = read_ints_from_line(current_line);
    }
}

void InputParser::parse_forces()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# FORCES") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## VOLUME FORCE") != 0 && current_line.compare("# END FORCES") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## VOLUME FORCE") == 0)
    {
        inputs.has_forces = true;
        inputs.has_volume_force = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.volume_force = read_matrix_row_from_line(current_line);
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("# END FORCES") != 0 && current_line.compare("# END FILE") != 0)
    {
        while (current_line.compare("## SURFACE FORCE") != 0 && current_line.compare("# END FORCES") != 0 && current_line.compare("# END FILE") != 0)
        {
            new_idx++;
            current_line = lines[new_idx];
        }
        if (current_line.compare("## SURFACE FORCE") == 0)
        {
            inputs.has_forces = true;
            inputs.has_surface_forces = true;
            new_idx++;
            current_line = lines[new_idx];
            inputs.tags_surface_forces.push_back(std::stoi(current_line));
            new_idx++;
            current_line = lines[new_idx];
            inputs.surface_forces.push_back(read_matrix_row_from_line(current_line));
        }
    }
}

void InputParser::parse_damping()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# DAMPING") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## MASS") != 0 && current_line.compare("# END DAMPING") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## MASS") == 0)
    {
        inputs.has_damping = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.damping_alpha_M = std::stof(current_line);
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("## STIFFNESS") != 0 && current_line.compare("# END DAMPING") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## STIFFNESS") == 0)
    {
        inputs.has_damping = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.damping_alpha_K = std::stof(current_line); 
    }
}

void InputParser::parse_solver()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# SOLVER") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    if (current_line.compare("# SOLVER") == 0)
    {
        inputs.has_solver = true;
        int new_idx = idx;
        new_idx++;
        current_line = lines[new_idx];
        inputs.solver_type = current_line;
        if (inputs.solver_type.compare("MODAL") == 0)
        {
            while (current_line.compare("## N MODES") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## N MODES") == 0)
            {
                new_idx++;
                current_line = lines[new_idx];
                inputs.modal_n_modes = std::stoi(current_line); 
            }
        }
        else if (inputs.solver_type.compare("FREQUENCYSWEEP") == 0)
        {
            while (current_line.compare("## COMPUTE ROM") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## COMPUTE ROM") == 0)
            {
                inputs.solver_computes_rom = true;
                new_idx++;
                current_line = lines[new_idx];
                inputs.frequency_n_modes = std::stoi(current_line); 
            }
            new_idx = idx;
            current_line = lines[new_idx];
            while (current_line.compare("## LOAD ROM") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## LOAD ROM") == 0)
            {
                inputs.solver_computes_rom = false;
                new_idx++;
                current_line = lines[new_idx];
                inputs.frequency_path_to_basis = current_line; 
            }
            new_idx = idx;
            current_line = lines[new_idx];
            while (current_line.compare("## MIN FREQUENCY") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## MIN FREQUENCY") == 0)
            {
                new_idx++;
                current_line = lines[new_idx];
                inputs.frequency_min_freq = std::stof(current_line); 
            }
            new_idx = idx;
            current_line = lines[new_idx];
            while (current_line.compare("## MAX FREQUENCY") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## MAX FREQUENCY") == 0)
            {
                new_idx++;
                current_line = lines[new_idx];
                inputs.frequency_max_freq = std::stof(current_line); 
            }
            new_idx = idx;
            current_line = lines[new_idx];
            while (current_line.compare("## N FREQUENCIES") != 0 && current_line.compare("# END SOLVER") != 0 && current_line.compare("# END FILE") != 0)
            {
                new_idx++;
                current_line = lines[new_idx];
            }
            if (current_line.compare("## N FREQUENCIES") == 0)
            {
                new_idx++;
                current_line = lines[new_idx];
                inputs.frequency_n_freq = std::stoi(current_line); 
            }
        }
    }
}

void InputParser::parse_output()
{
    int idx = 0;
    std::string current_line = lines[idx];
    while (current_line.compare("# OUTPUT") != 0 && current_line.compare("# END FILE") != 0)
    {
        idx++;
        current_line = lines[idx];
    }
    int new_idx = idx;
    while (current_line.compare("## STUDY NAME") != 0 && current_line.compare("# END OUTPUT") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## STUDY NAME") == 0)
    {
        new_idx++;
        current_line = lines[new_idx];
        inputs.output_name = current_line;
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("## PATH") != 0 && current_line.compare("# END OUTPUT") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## PATH") == 0)
    {
        inputs.has_output = true;
        new_idx++;
        current_line = lines[new_idx];
        inputs.output_path = current_line;
    }
    new_idx = idx;
    current_line = lines[new_idx];
    while (current_line.compare("## PLOTTED MODE NUM") != 0 && current_line.compare("# END OUTPUT") != 0 && current_line.compare("# END FILE") != 0)
    {
        new_idx++;
        current_line = lines[new_idx];
    }
    if (current_line.compare("## PLOTTED MODE NUM") == 0)
    {
        new_idx++;
        current_line = lines[new_idx];
        inputs.plotted_mode_num = std::stoi(current_line); 
    }
}

void InputParser::parse_input_file()
{
    inputs.has_mesh = false;
    inputs.has_material = false;
    inputs.has_dirichlet = false;
    inputs.has_forces = false;
    inputs.has_surface_forces = false;
    inputs.has_volume_force = false;
    inputs.has_damping = false;
    inputs.has_solver = false;
    inputs.solver_computes_rom = false;
    inputs.has_output = false;

    std::ifstream input_file(path_to_file);
    std::string current_line;

    current_line = "init";
    get_lines(input_file, current_line);

    parse_mesh();
    parse_material();
    parse_dirichlet();
    parse_forces();
    parse_damping();
    parse_solver();
    parse_output();

    input_file.close();
}