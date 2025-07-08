//
// mesh.h
//

#ifndef mesh_h
#define mesh_h

#include <vector>

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

        // Imports mesh data (node coordinates and connectivity tables) from a GMSH tetrahedral mesh file saved in Matlab format 
        // (uncheck the "Save all elements" box in order to preserve physical group tags).
        // Only elements affected to a physical group of dimension 2 (triangles) and 3 (tetrahedra) will be read.
        void import_gmsh_matlab();
};

#endif