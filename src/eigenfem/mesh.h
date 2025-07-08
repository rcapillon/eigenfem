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

struct Element {
    int tag;
    std::vector<int> nodes_tags;

    Element(int t, std::vector<int> nt) : tag(t), nodes_tags(nt) {}
};

struct ElementsTable {
    public:
        int tag;
        int dim;
        std::vector<Element> elements;
        int n_elements = elements.size();
        Eigen::MatrixXf table;

    ElementsTable(int t, int d, std::vector<Element> e) : tag(t), dim(d), elements(e) {
        table = Eigen::MatrixXf(n_elements, dim);
    }
};

class Mesh {
    public:
        Mesh(std::string path_to_mesh);
        ~Mesh();

        std::string mesh_path;
        std::vector<Node> vec_Nodes;
        std::vector<ElementsTable> tables_Elements;
        int num_3D_elements_table;

        void import_gmsh();
};

#endif