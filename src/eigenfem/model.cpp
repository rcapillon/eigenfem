//
// model.cpp
//

#include "model.h"


Model::Model(Mesh msh,
    std::vector<int> dir_tags, 
    std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces,
    std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces,
    float a_M, float a_K) : alpha_M{0.}, alpha_K{0.}
{
    mesh = msh;
    dirichlet_tags = dir_tags;
    surface_forces = surf_forces;
    volume_forces = vol_forces;
    alpha_M = a_M;
    alpha_K = a_K;
};

void Model::apply_dirichlet()
{
    
}