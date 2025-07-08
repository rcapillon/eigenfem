//
// mesh.cpp
//

#include "mesh.h"


void Mesh::import_gmsh()
{
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
	Mesh::nodes_table = NodesTable(vec_Nodes.size());
	for (size_t k = 0; k < vec_Nodes.size(); k++)
	{
		Mesh::nodes_table.table(k, 0) = vec_Nodes[k].coords[0];
		Mesh::nodes_table.table(k, 1) = vec_Nodes[k].coords[1];
		Mesh::nodes_table.table(k, 2) = vec_Nodes[k].coords[2];
		Mesh::nodes_table.table(k, 3) = vec_Nodes[k].tag;
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
		if (table_dim == 3)
		{
			Mesh::num_3D_elements_table = int(i);
		}
		

		int n_elements = elements.entity_blocks[i].data.size();
		std::vector<Element> vec_elements;
		for (size_t j = 0; j < n_elements; j++)
		{
			int element_tag = elements.entity_blocks[i].data[(table_dim + 2) * i];
			std::vector<int> nodes_tags;
			for (size_t k = 1; k < table_dim + 2; k++)
			{
				nodes_tags.push_back(elements.entity_blocks[i].data[(table_dim + 2) * i + k]);
			}
			Element element(element_tag, nodes_tags);	
			vec_elements.push_back(element);		
		}
		ElementsTable tmp_elements_table(table_tag, table_dim, vec_elements);

		for (size_t m = 0; m < vec_elements.size(); m++)
		{
			Element element = vec_elements[m];
			for (size_t n = 0; n < element.nodes_tags.size(); n++)
			{
				int node_tag = element.nodes_tags[n];
				bool found = false;
				int count = 0;
				while (!found)
				{
					int nodes_node_tag = Mesh::nodes_table.table(count, 3);
					if (node_tag == nodes_node_tag)
					{
						found = true;
					}
					else
					{
						count++;
					}
				}
				tmp_elements_table.table(m, n) = count;
			}
			
		}
		Mesh::elements_tables.push_back(tmp_elements_table);
	}
}
