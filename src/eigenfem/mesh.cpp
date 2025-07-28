//
// mesh.cpp
//

#include "mesh.h"


// Reads lines in Matlab-format GMSH mesh files for node numbers in connectivity.
// Ignores physical group tag and adapts the numbering to start at 0.
std::vector<int> read_node_nums_from_line(std::string line)
{
	std::vector<int> nums;
	size_t prev = 0;
	size_t pos;
	std::string substr;
    while ((pos = line.find_first_of(" ;", prev)) != std::string::npos)
    {
        if (pos > prev)
            nums.push_back(std::stoi(line.substr(prev, pos-prev)));
        prev = pos + 1;
    }
    if (prev < line.length())
        nums.push_back(std::stoi(line.substr(prev, std::string::npos)));

	for (size_t i = 0; i < nums.size() - 1; i++)
	{
		nums[i] -= 1;
	}
	

	return nums;
}

// Reads lines in Matlab-format GMSH mesh files for node coordinates in 3D space.
std::vector<float> read_node_coords_from_line(std::string line)
{
	std::vector<float> coords;
	size_t prev = 0;
	size_t pos;
	std::string substr;
    while ((pos = line.find_first_of(" ;", prev)) != std::string::npos)
    {
        if (pos > prev)
            coords.push_back(std::stof(line.substr(prev, pos-prev)));
        prev = pos + 1;
    }
    if (prev < line.length())
        coords.push_back(std::stof(line.substr(prev, std::string::npos)));

	return coords;
}

// Reads lines in Matlab-format GMSH mesh files for nodes groups (0D physical groups)
std::vector<int> read_nodes_groups_line(std::string line)
{
	std::vector<int> node_num_and_tag;
	size_t prev = 0;
	size_t pos;
	std::string substr;
	while ((pos = line.find_first_of(" ", prev)) != std::string::npos)
    {
        if (pos > prev)
            node_num_and_tag.push_back(std::stoi(line.substr(prev, pos-prev)));
        prev = pos + 1;
    }
    if (prev < line.length())
        node_num_and_tag.push_back(std::stoi(line.substr(prev, std::string::npos)));

	return node_num_and_tag;
}

// Mesh class constructors
Mesh::Mesh(std::string path_to_mesh) 
{
	mesh_path = path_to_mesh;
	material = Material(0., 0., 0.);
};

Mesh::Mesh(std::string path_to_mesh, Material mat) 
{
	mesh_path = path_to_mesh;
	material = mat;
};

// Imports mesh data (node coordinates and connectivity tables) from a GMSH tetrahedral (1st order) mesh file 
// saved in Matlab format.
// Uncheck the "Save all elements" checkbox in order to preserve physical group tags.
// Only elements affected to a physical group of dimension 2 (triangles) and 3 (tetrahedra) will be read.
void Mesh::import_gmsh_matlab()
{
	// open file
	std::ifstream file(mesh_path);
    std::string current_line; 

	// initialize file reading and get number of nodes in the mesh
	for (size_t i = 0; i < 5; i++)
	{
		std::getline(file, current_line);
	}
	n_nodes = std::stoi(current_line.substr(12, std::string::npos));
	n_dofs = 3 * n_nodes;
	std::getline(file, current_line);
	std::getline(file, current_line);

	// fill table of node coordinates
	table_nodes = Eigen::MatrixXf(n_nodes, 3);
	int n_row = 0;
	while (current_line.compare("];") != 0)
	{
		std::vector<float> coords = read_node_coords_from_line(current_line);
		table_nodes(n_row, 0) = coords[0];
		table_nodes(n_row, 1) = coords[1];
		table_nodes(n_row, 2) = coords[2];
		n_row++;
		std::getline(file, current_line);
	}

	// initialize triangle connectivity reading
	for (size_t i = 0; i < 4; i++)
	{
		std::getline(file, current_line);
	}

	// fill tables of triangle connectivity by physical group tag
	std::vector<int> nums = read_node_nums_from_line(current_line);
	int previous_tag = nums[nums.size() - 1];
	int current_tag = previous_tag;
	std::vector<int> current_nums = std::vector<int>(nums.begin(), nums.end() - 1);
	std::vector<std::vector<int>> current_tri_group;
	current_tri_group.push_back(current_nums);

	std::getline(file, current_line);
	while (current_line.compare("];") != 0)
	{
		std::vector<int> nums = read_node_nums_from_line(current_line);
		current_tag = nums[nums.size() - 1];

		if (std::find(tris_tags.begin(), tris_tags.end(), current_tag) == tris_tags.end())
		{
			tris_tags.push_back(current_tag);
		}
		
		while (current_tag == previous_tag)
		{
			std::vector<int> current_nums = std::vector<int>(nums.begin(), nums.end() - 1);
			current_tri_group.push_back(current_nums);

			std::getline(file, current_line);
			if (current_line.compare("];") != 0)
			{
				nums = read_node_nums_from_line(current_line);
				previous_tag = current_tag;
				current_tag = nums[nums.size() - 1];
			}
			else
			{
				current_tag = previous_tag + 1;
			}
		}
		
		int group_n_elements = current_tri_group.size();
		Eigen::MatrixXi current_table_tri(group_n_elements, 3);
		for (size_t i = 0; i < group_n_elements; i++)
		{
			current_table_tri(i, 0) = current_tri_group[i][0];
			current_table_tri(i, 1) = current_tri_group[i][1];
			current_table_tri(i, 2) = current_tri_group[i][2];
		}
		tables_tris.push_back(current_table_tri);
		current_tri_group.clear();

		previous_tag = current_tag;
	}

	// initialize tetrahedra reading
	std::getline(file, current_line);
	std::getline(file, current_line);

	// fill table of tetrahedra connectivity and create elements
	std::vector<std::vector<int>> nums_tetra;
	while (current_line.compare("];") != 0)
	{
		std::vector<int> nums = read_node_nums_from_line(current_line);
		nums_tetra.push_back(nums);
		std::getline(file, current_line);
	}
	n_elements = nums_tetra.size();
	table_tets = Eigen::MatrixXi(n_elements, 4);
	for (size_t i = 0; i < n_elements; i++)
	{
		table_tets(i, 0) = nums_tetra[i][0];
		table_tets(i, 1) = nums_tetra[i][1];
		table_tets(i, 2) = nums_tetra[i][2];
		table_tets(i, 3) = nums_tetra[i][3];

		Eigen::MatrixXf nodes_coords(4, 3);
		std::vector<int> inds = {0, 1, 2};
		nodes_coords(0, inds) = table_nodes(nums_tetra[i][0], inds);
		nodes_coords(1, inds) = table_nodes(nums_tetra[i][1], inds);
		nodes_coords(2, inds) = table_nodes(nums_tetra[i][2], inds);
		nodes_coords(3, inds) = table_nodes(nums_tetra[i][3], inds);
		Element element = Element(int(i), material, nums_tetra[i], nodes_coords);
		elements.push_back(element);
	}

	// initialize nodes groups reading
	std::getline(file, current_line);
	std::getline(file, current_line);

	// fill nodes groups and save tags
	if (!current_line.empty())
	{
		std::vector<int> all_nodes_in_groups;
		std::vector<int> all_nodes_tags;
		while (current_line.compare("];") != 0)
		{
			std::vector<int> node_num_and_tag = read_nodes_groups_line(current_line);
			all_nodes_in_groups.push_back(node_num_and_tag[0]);
			all_nodes_tags.push_back(node_num_and_tag[1]);
			std::getline(file, current_line);
		}
		std::vector<int> all_nodes_tags_copy = all_nodes_tags ;
		int n_nodes_tags = std::unique(all_nodes_tags_copy.begin(), all_nodes_tags_copy.end()) - all_nodes_tags_copy.begin();
		all_nodes_tags_copy.erase(std::unique(all_nodes_tags_copy.begin(), all_nodes_tags_copy.end()), all_nodes_tags_copy.end());
		for (size_t i = 0; i < n_nodes_tags; i++)
		{
			std::vector<int> nodes_in_group;
			nodes_groups.push_back(nodes_in_group);
		}
		for (size_t i = 0; i < n_nodes_tags; i++)
		{
			int current_tag = all_nodes_tags_copy[i];
			nodes_tags.push_back(current_tag);
			for (size_t j = 0; j < all_nodes_in_groups.size(); j++)
			{
				int tested_tag = all_nodes_tags[j];
				if (tested_tag == current_tag)
				{
					nodes_groups[i].push_back(all_nodes_in_groups[j]);
				}	
			}
		}
	}
	
	file.close();
}
