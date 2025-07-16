//
// model.cpp
//

#include "model.h"


Model::Model(Mesh msh, std::vector<int> dir_tags)
{
    mesh = msh;
    dirichlet_tags = dir_tags;
    alpha_M = 0.;
    alpha_K = 0.;
};

Model::Model(Mesh msh, std::vector<int> dir_tags, float a_M, float a_K)
{
    mesh = msh;
    dirichlet_tags = dir_tags;
    alpha_M = a_M;
    alpha_K = a_K;
};

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

void Model::create_dof_lists()
{
    for (size_t i = 0; i < dirichlet_tags.size(); i++)
    {
        auto it = std::find(mesh.tris_tags.begin(), mesh.tris_tags.end(), dirichlet_tags[i]);
        int idx_tables_tris = int(it - mesh.tris_tags.begin());
        for (size_t j = 0; j < mesh.tables_tris[idx_tables_tris].rows(); j++)
        {
            for (size_t k = 0; k < 3; k++)
            {
                int dof1_num = mesh.tables_tris[idx_tables_tris](j, k) * 3;
                if (std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), dof1_num) != dirichlet_dofs.end())
                {
                dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3);
                dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3 + 1);
                dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3 + 2);
                }
            }
        }
    }
    
};

void Model::apply_dirichlet() {};