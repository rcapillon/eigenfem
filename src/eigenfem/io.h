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
std::vector<int> read_ints_from_line(std::string line);

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

struct Inputs
{
    std::string mesh_path;

    float material_rho;
    float material_youngmodulus;
    float material_poissonratio;

    std::vector<int> dirichlet_tags;

    std::vector<float> volume_force;
    std::vector<int> tags_surface_forces;
    std::vector<std::vector<float>> surface_forces;

    float damping_alpha_M;
    float damping_alpha_K;

    std::string solver_type;

    int modal_n_modes;

    int frequency_n_modes;
    std::string frequency_path_to_basis;
    float frequency_min_freq;
    float frequency_max_freq;
    int frequency_n_freq;

    std::string output_name;
    std::string output_path;

    bool has_surface_force;
    bool has_volume_force;
    bool has_damping;
    bool compute_rom;
    bool has_output;

    Inputs() {};
};

class InputParser
{
    public:
        InputParser() {};
        InputParser(std::string path_to_input_file);
        ~InputParser() {};

        std::string path_to_file;

        std::string current_line;

        Inputs inputs;

        void parse_mesh(std::ifstream& file);
        void parse_material(std::ifstream& file);
        void parse_dirichlet(std::ifstream& file);
        void parse_forces(std::ifstream& file);
        void parse_damping(std::ifstream& file);
        void parse_solver(std::ifstream& file);
        void parse_output(std::ifstream& file);
        void parse_input_file();
};

#endif