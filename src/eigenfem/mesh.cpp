//
// mesh.cpp
//

#include <iostream>
#include <string>

#include "../../third-party/MshIO-0.0.1/include/mshio/mshio.h"


struct Node
{
	int tag;
	std::vector<float> coords;

	Node(int t, std::vector<float> c) : tag(t), coords(c) {}
};

class Mesh
{
public:
	Mesh();
	~Mesh();

	string mesh_path;
	std::vector<Node> vec_Nodes;

	
};