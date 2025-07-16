//
// model.cpp
//

#include <iostream>

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
            for (size_t k = 0; k < mesh.tables_tris[idx_tables_tris].cols(); k++)
            {
                int dof1_num = mesh.tables_tris[idx_tables_tris](j, k) * 3;
                if (std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), dof1_num) == dirichlet_dofs.end())
                {
                    dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3);
                    dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3 + 1);
                    dirichlet_dofs.push_back(mesh.tables_tris[idx_tables_tris](j, k) * 3 + 2);
                }
            }
        }
    }

    for (size_t i = 0; i < mesh.table_tets.rows(); i++)
    {
        for (size_t j = 0; j < mesh.table_tets.cols(); j++)
        {
            int dof1_num = mesh.table_tets(i, j) * 3;
            if (std::find(free_dofs.begin(), free_dofs.end(), dof1_num) == free_dofs.end() 
            && std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), dof1_num) == dirichlet_dofs.end())
            {
                free_dofs.push_back(mesh.table_tets(i, j) * 3);
                free_dofs.push_back(mesh.table_tets(i, j) * 3 + 1);
                free_dofs.push_back(mesh.table_tets(i, j) * 3 + 2);
            }
        }
    }
};

void Model::assemble_M()
{
    SpMat mat_M(mesh.n_dofs, mesh.n_dofs);
    std::vector<triplet> triplets;

    for (size_t i = 0; i < mesh.n_elements; i++)
    {
        std::vector<int> ind_I;
        std::vector<int> ind_J;
        for (size_t j = 0; j < mesh.elements[i].n_dofs; j++)
        {
            for (size_t k = 0; k < mesh.elements[i].n_dofs; k++)
            {
                ind_I.push_back(k);
                ind_J.push_back(j);
            }
            
        }
        Eigen::VectorXf flattened_mat_Me = mesh.elements[i].compute_mat_Me().reshaped();

        for (size_t m = 0; m < ind_I.size(); m++)
        {
            triplets.push_back(triplet(ind_I[m], ind_J[m], flattened_mat_Me(m)));
        }
    }
    mat_M.setFromTriplets(triplets.begin(), triplets.end());
};

void Model::assemble_K()
{
    SpMat mat_K(mesh.n_dofs, mesh.n_dofs);
    std::vector<triplet> triplets;

    for (size_t i = 0; i < mesh.n_elements; i++)
    {
        std::vector<int> ind_I;
        std::vector<int> ind_J;
        for (size_t j = 0; j < mesh.elements[i].n_dofs; j++)
        {
            for (size_t k = 0; k < mesh.elements[i].n_dofs; k++)
            {
                ind_I.push_back(k);
                ind_J.push_back(j);
            }
            
        }
        Eigen::VectorXf flattened_mat_Ke = mesh.elements[i].compute_mat_Ke().reshaped();

        for (size_t m = 0; m < ind_I.size(); m++)
        {
            triplets.push_back(triplet(ind_I[m], ind_J[m], flattened_mat_Ke(m)));
        }
    }
    mat_K.setFromTriplets(triplets.begin(), triplets.end());
};

void Model::assemble_M_K()
{
    SpMat mat_M(mesh.n_dofs, mesh.n_dofs);
    SpMat mat_K(mesh.n_dofs, mesh.n_dofs);
    std::vector<triplet> triplets_M;
    std::vector<triplet> triplets_K;

    for (size_t i = 0; i < mesh.n_elements; i++)
    {
        std::vector<int> ind_I;
        std::vector<int> ind_J;
        for (size_t j = 0; j < mesh.elements[i].n_dofs; j++)
        {
            for (size_t k = 0; k < mesh.elements[i].n_dofs; k++)
            {
                ind_I.push_back(k);
                ind_J.push_back(j);
            }
            
        }
        Eigen::VectorXf flattened_mat_Me = mesh.elements[i].compute_mat_Me().reshaped();
        Eigen::VectorXf flattened_mat_Ke = mesh.elements[i].compute_mat_Ke().reshaped();

        for (size_t m = 0; m < ind_I.size(); m++)
        {
            triplets_M.push_back(triplet(ind_I[m], ind_J[m], flattened_mat_Me(m)));
            triplets_K.push_back(triplet(ind_I[m], ind_J[m], flattened_mat_Ke(m)));
        }
    }
    mat_M.setFromTriplets(triplets_M.begin(), triplets_M.end());
    mat_K.setFromTriplets(triplets_K.begin(), triplets_K.end());
};

void Model::compute_D_Rayleigh()
{
    mat_D = alpha_M * mat_M + alpha_K * mat_K;
};

void Model::apply_dirichlet()
{
    if (mat_M.rows() * mat_M.cols() != 0)
    {
        mat_Mff = double_slice_spmat(mat_M, free_dofs, free_dofs);
    }
    if (mat_K.rows() * mat_K.cols() != 0)
    {
        mat_Kff = double_slice_spmat(mat_K, free_dofs, free_dofs);
    }
    if (mat_D.rows() * mat_D.cols() != 0)
    {
        mat_Dff = double_slice_spmat(mat_D, free_dofs, free_dofs);
    }
};