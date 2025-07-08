//
// main.cpp
//

#include <iostream>

#include "../../third-party/MshIO-0.0.1/include/mshio/mshio.h"
#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"
#include "../../third-party/spectra-1.1.0/include/Spectra/GenEigsSolver.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseGenMatProd.h"

#include "mesh.h"


int main()
{
    mshio::MshSpec spec = mshio::load_msh("../data/bar.msh");

    auto& elements = spec.elements;

    std::cout << "Number of element blocks:\n" << elements.entity_blocks.size() << "\n";

    std::cout << "Tag and dimension of each element block:\n";
    for (int i = 0; i < elements.entity_blocks.size(); i++)
    {
        std::cout << "Tag " << elements.entity_blocks[i].entity_tag << " , dimension " << elements.entity_blocks[i].entity_dim << "\n";
    }

    std::cout << "Size of element data of each element block:\n";
    for (int i = 0; i < elements.entity_blocks.size(); i++)
    {
        std::cout << "Tag " << elements.entity_blocks[i].entity_tag << " , size of element data: " << elements.entity_blocks[i].data.size() << "\n";
    }

    auto& nodes = spec.nodes;

    std::cout << "Number of node blocks:\n" << nodes.entity_blocks.size() << "\n";

    int total_nodes = 0;
    std::cout << "Dimension and number of nodes for each node block:\n";
    for (int i = 0; i < nodes.entity_blocks.size(); i++)
    {
        std::cout << "Tag " << nodes.entity_blocks[i].entity_tag << " , dimension " << nodes.entity_blocks[i].entity_dim << " , number of nodes: " << nodes.entity_blocks[i].num_nodes_in_block << "\n";
        total_nodes += nodes.entity_blocks[i].num_nodes_in_block;
    }

    std::cout << "Total number of nodes: " << total_nodes << "\n";

    std::cout << "Size of data in last node block:\n";
    std::cout << nodes.entity_blocks[26].data.size() << "\n";

    std::cout << "Coordinates of first node in last node block:\n";
    for (int i = 0; i < 3; i++)
    {
        std::cout << nodes.entity_blocks[26].data[i] << " " ;
    }
    std::cout << "\n";

    return 0;
}