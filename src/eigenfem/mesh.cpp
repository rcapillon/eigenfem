//
// mesh.cpp
//

#include <fstream>
#include <string>
#include <iostream>

#include "mesh.h"


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

void Mesh::import_gmsh_matlab()
{
	// open file
	std::ifstream file(Mesh::mesh_path);
    std::string current_line; 

	// initialize file reading and get number of nodes in the mesh
	for (size_t i = 0; i < 5; i++)
	{
		std::getline(file, current_line);
	}
	Mesh::n_nodes = std::stoi(current_line.substr(12, std::string::npos));
	std::getline(file, current_line);
	std::getline(file, current_line);

	// fill table of node coordinates
	Mesh::table_nodes = Eigen::MatrixXf(Mesh::n_nodes, 3);
	int n_row = 0;
	while (current_line.compare("];") != 0)
	{
		std::vector<float> coords = read_node_coords_from_line(current_line);
		Mesh::table_nodes(n_row, 0) = coords[0];
		Mesh::table_nodes(n_row, 1) = coords[1];
		Mesh::table_nodes(n_row, 2) = coords[2];
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
		Eigen::MatrixXf current_table_tri(group_n_elements, 3);
		for (size_t i = 0; i < group_n_elements; i++)
		{
			current_table_tri(i, 0) = current_tri_group[i][0];
			current_table_tri(i, 1) = current_tri_group[i][1];
			current_table_tri(i, 2) = current_tri_group[i][2];
		}
		Mesh::tables_tris.push_back(current_table_tri);

		previous_tag = current_tag;
	}

	file.close();
}
