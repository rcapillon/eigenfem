//
// io.h
//

#ifndef io_h
#define io_h

#include <fstream>
#include <string>
#include <algorithm>

#include "../../third-party/eigen-3.4.0/Eigen/Core"

#include "mesh.h"


std::vector<float> read_matrix_row_from_line(std::string line);

class VTKwriter
{
    public:
        VTKwriter() {};
        VTKwriter(Mesh msh, Eigen::VectorXf vec_displacement);
        VTKwriter(Mesh msh, Eigen::MatrixXf mat_displacement);
        ~VTKwriter() {};

        Mesh mesh;
        Eigen::VectorXf U;
        Eigen::MatrixXf mat_U;

        Mesh deformed_mesh;

        void add_U_to_mesh();
        void add_U_to_mesh(int idx_U);
        void reset_deformed_mesh();
        void write_deformed_mesh(std::string path_to_folder, std::string filename_wo_extension, std::string trail = "");
        void write_mesh_animation(std::string path_to_folder, std::string filename_wo_extension);
};

class DATio
{
    public:
        DATio() {};
        DATio(Eigen::MatrixXf mat);
        ~DATio() {};

        Eigen::MatrixXf matrix;

        void write_matrix(std::string path_to_file);
        Eigen::MatrixXf load_dat(std::string path_to_file);
};

class InputParser
{
    public:
        InputParser() {};
        InputParser(std::string path_to_input_file);
        ~InputParser() {};
};

#endif