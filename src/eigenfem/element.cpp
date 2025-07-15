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

void Element::compute_jacobian_at_gauss_points() {};

Eigen::MatrixXf Element::compute_mat_Me() {};
Eigen::MatrixXf Element::compute_mat_Ke() {};