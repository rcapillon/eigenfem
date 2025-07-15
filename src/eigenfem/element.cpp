//
// element.cpp
//

#include "element.h"


Eigen::MatrixXf mat_G = compute_mat_G();
Eigen::MatrixXf mat_P = compute_mat_P();
Eigen::MatrixXf precalc_gauss_points = get_gauss_points();
Eigen::MatrixXf precalc_gauss_weights = get_gauss_weights();
Eigen::MatrixXf precalc_nodes_reference_coords = get_reference_coords();
Eigen::MatrixXf precalc_shapefun_coeffs = get_shapefun_coeffs();

MatsGauss precalc_mats_gauss = compute_mats_gauss();

Element::Element(int num, Material mat, std::vector<int> nodes_numbers, Eigen::MatrixXf nodes_coordinates)
{
    number = num;
    material = mat;
    nodes_num = nodes_numbers;
    nodes_coords = nodes_coordinates;

    n_nodes = 4;
    n_dofs = 3 * n_nodes;
    nodes_reference_coords = precalc_nodes_reference_coords;
    gauss_points = precalc_gauss_points;
    gauss_weights = precalc_gauss_weights;
    for (size_t i = 0; i < nodes_num.size(); i++)
    {
        dofs_num.push_back(3 * nodes_num[i]);
        dofs_num.push_back(3 * nodes_num[i] + 1);
        dofs_num.push_back(3 * nodes_num[i] + 2);
    }
    for (size_t i = 0; i < nodes_coords.rows(); i++)
    {
        vec_nodes_coords[3 * i] = nodes_coords[i, 0];
        vec_nodes_coords[3 * i + 1] = nodes_coords[i, 1];
        vec_nodes_coords[3 * i + 2] = nodes_coords[i, 2];
    }
    mats_gauss = precalc_mats_gauss;
};

MatsJ Element::compute_jacobian_at_gauss_points(int gauss_idx)
{
    std::vector<int> inds0 = {0, 1, 2};
    std::vector<int> inds1 = {3, 4, 5};
    std::vector<int> inds2 = {6, 7, 8};

    Eigen::MatrixXf mat_Ddx = mats_gauss.vec_mat_De_gauss[gauss_idx](inds0, Eigen::placeholders::all);
    Eigen::MatrixXf mat_Ddy = mats_gauss.vec_mat_De_gauss[gauss_idx](inds1, Eigen::placeholders::all);
    Eigen::MatrixXf mat_Ddz = mats_gauss.vec_mat_De_gauss[gauss_idx](inds2, Eigen::placeholders::all);

    Eigen::MatrixXf vec_J1 = mat_Ddx * vec_nodes_coords;
    Eigen::MatrixXf vec_J2 = mat_Ddy * vec_nodes_coords;
    Eigen::MatrixXf vec_J3 = mat_Ddz * vec_nodes_coords;

    Eigen::MatrixXf mat_J(3, 3);
    mat_J(0, inds0) = vec_J1.transpose();
    mat_J(1, inds0) = vec_J2.transpose();
    mat_J(2, inds0) = vec_J3.transpose();

    float det_J = mat_J.determinant();

    Eigen::MatrixXf mat_I = Eigen::MatrixXf::Identity(3, 3);
    Eigen::MatrixXf mat_invJ = mat_J.ldlt().solve(mat_I);
    Eigen::MatrixXf mat_invJJJ(9, 9);
    mat_invJJJ(inds0, inds0) = mat_invJ;
    mat_invJJJ(inds1, inds1) = mat_invJ;
    mat_invJJJ(inds2, inds2) = mat_invJ;

    MatsJ mats_J(det_J, mat_invJJJ);
};

Eigen::MatrixXf Element::compute_mat_Me()
{
    Eigen::MatrixXf mat_Me(n_dofs, n_dofs);
    for (size_t i = 0; i < gauss_points.rows(); i++)
    {
        MatsJ mats_J = compute_jacobian_at_gauss_points(i);
        mat_Me += gauss_weights(i) * material.rho * mats_J.det_J * mats_gauss.vec_mat_EeTEe_gauss[i];
    }

    return mat_Me;
};

Eigen::MatrixXf Element::compute_mat_Ke()
{
    Eigen::MatrixXf mat_Ke(n_dofs, n_dofs);
    for (size_t i = 0; i < gauss_points.rows(); i++)
    {
        MatsJ mats_J = compute_jacobian_at_gauss_points(i);
        Eigen::MatrixXf mat_B = mat_G * mats_J.mat_invJJJ * mat_P * mats_gauss.vec_mat_De_gauss[i];
        mat_Ke += gauss_weights(i) * mats_J.det_J * mat_B.transpose() * material.mat_C * mat_B;
    }
    
};
