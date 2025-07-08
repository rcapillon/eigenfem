//
// mesh.cpp
//

#include "mesh.h"


Mesh::Mesh(std::string path_to_mesh) {
	mesh_path = path_to_mesh;
};

void Mesh::import_gmsh() {
	//
	// import all nodes
	//
	mshio::MshSpec spec = mshio::load_msh(mesh_path);
	auto& nodes = spec.nodes;
	size_t n_nodeblocks = nodes.entity_blocks.size();
	for (size_t i = 0; i < n_nodeblocks; i++)
	{
		size_t n_nodes_in_block = nodes.entity_blocks[i].num_nodes_in_block;
		for (size_t j = 0; j < n_nodeblocks; j++)
		{
			int node_tag = nodes.entity_blocks[i].tags[j];
			float node_x = nodes.entity_blocks[i].data[3 * j];
			float node_y = nodes.entity_blocks[i].data[3 * j + 1];
			float node_z = nodes.entity_blocks[i].data[3 * j + 2];
			std::vector<float> node_coords = {node_x, node_y, node_z};

			Node node(node_tag, node_coords);
			vec_Nodes.push_back(node);
		}
	}

	//
	// import all elements
	//
	auto& elements = spec.elements;
	size_t n_elementblocks = elements.entity_blocks.size();
	for (size_t i = 0; i < n_elementblocks; i++)
	{
		int table_tag = elements.entity_blocks[i].entity_tag;
		int table_dim = elements.entity_blocks[i].entity_dim;
	}
	

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
}
