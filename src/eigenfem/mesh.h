//
// mesh.h
//

#ifndef mesh_h
#define mesh_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"


class Mesh {
    public:
        Mesh(std::string path_to_mesh) {
            Mesh::mesh_path = path_to_mesh;
        };
        ~Mesh() {};

        std::string mesh_path;
        int n_nodes;
        int n_elements;
        Eigen::MatrixXf table_nodes;
        Eigen::MatrixXf table_tets;
        std::vector<Eigen::MatrixXf> tables_tris;

        void import_gmsh_matlab();
};

#endif