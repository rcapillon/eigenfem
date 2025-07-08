//
// mesh.h
//

#ifndef mesh_h
#define mesh_h

#include <string>
#include <vector>

#include "../../third-party/MshIO-0.0.1/include/mshio/mshio.h"
#include "../../third-party/eigen-3.4.0/Eigen/Core"

struct Node {
    int tag;
	std::vector<float> coords;

    Node(int t, std::vector<float> c) : tag(t), coords(c) {}
};

struct NodesTable {
    int n_nodes;
    Eigen::MatrixXf table;

    NodesTable() {};
    NodesTable(int n) : n_nodes(n) {
        table = Eigen::MatrixXf(n_nodes, 4);
    }
};

struct Element {
    int tag;
    std::vector<int> nodes_tags;

    Element(int t, std::vector<int> nt) : tag(t), nodes_tags(nt) {}
};

struct ElementsTable {
    int tag;
    int dim;
    std::vector<Element> elements;
    int n_elements = elements.size();
    Eigen::MatrixXf table;

    ElementsTable() {};
    ElementsTable(int t, int d, std::vector<Element> e) : tag(t), dim(d), elements(e) {
        table = Eigen::MatrixXf(n_elements, dim + 1);
    }
};

class Mesh {
    public:
        Mesh(std::string path_to_mesh) {
            Mesh::mesh_path = path_to_mesh;
        };
        ~Mesh() {};

        std::string mesh_path;
        std::vector<Node> vec_Nodes;
        NodesTable nodes_table;
        std::vector<ElementsTable> elements_tables;
        int num_3D_elements_table;

        void import_gmsh();
};

#endif