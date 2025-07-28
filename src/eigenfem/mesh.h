//
// mesh.h
//

#ifndef mesh_h
#define mesh_h

#include <fstream>
#include <string>

#include "../../third-party/eigen-3.4.0/Eigen/Core"

#include "element.h"
#include "materials.h"


std::vector<int> read_node_nums_from_line(std::string line);
std::vector<float> read_node_coords_from_line(std::string line);
std::vector<int> read_nodes_groups_line(std::string line);

class Mesh 
{
    public:
        Mesh() {};
        Mesh(std::string path_to_mesh);
        Mesh(std::string path_to_mesh, Material mat);
        ~Mesh() {};

        std::string mesh_path;
        Material material;
        int n_nodes;
        int n_dofs;
        int n_elements;
        Eigen::MatrixXf table_nodes;
        Eigen::MatrixXi table_tets;
        std::vector<Eigen::MatrixXi> tables_tris;
        std::vector<int> tris_tags;
        std::vector<std::vector<int>> nodes_groups;
        std::vector<int> nodes_tags;
        std::vector<Element> elements;

        void import_gmsh_matlab();
};

#endif