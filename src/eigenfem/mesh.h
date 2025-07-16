//
// mesh.h
//

#ifndef mesh_h
#define mesh_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"

#include "element.h"
#include "materials.h"


std::vector<int> read_node_nums_from_line(std::string line);
std::vector<float> read_node_coords_from_line(std::string line);

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
        Eigen::MatrixXf table_tets;
        std::vector<Eigen::MatrixXf> tables_tris;
        std::vector<int> tris_tags;
        std::vector<Element> elements;

        void import_gmsh_matlab();
};

#endif