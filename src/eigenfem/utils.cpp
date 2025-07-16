//
// utils.cpp
//

#include "utils.h"


Eigen::MatrixXf compute_mat_G()
{
    float s = sqrt(2.) / 2.;
    Eigen::MatrixXf mat_G = Eigen::MatrixXf::Zero(6, 9);
    mat_G(0, 0) = 1.;
    mat_G(1, 4) = 1.;
    mat_G(2, 8) = 1.;
    mat_G(3, 5) = s;
    mat_G(3, 7) = s;
    mat_G(4, 2) = s;
    mat_G(4, 6) = s;
    mat_G(5, 1) = s;
    mat_G(5, 3) = s;

    return mat_G;
}

Eigen::MatrixXf compute_mat_P()
{
    Eigen::MatrixXf mat_P = Eigen::MatrixXf::Zero(9, 9);
    mat_P(0, 0) = 1.;
    mat_P(1, 3) = 1.;
    mat_P(2, 6) = 1.;
    mat_P(3, 1) = 1.;
    mat_P(4, 4) = 1.;
    mat_P(5, 7) = 1.;
    mat_P(6, 2) = 1.;
    mat_P(7, 5) = 1.;
    mat_P(8, 8) = 1.;

    return mat_P;
}

Eigen::MatrixXf get_reference_coords()
{
    Eigen::MatrixXf nodes_reference_coords = Eigen::MatrixXf::Zero(4, 3);
    nodes_reference_coords(0, 0) = 0.;
    nodes_reference_coords(0, 1) = 0.;
    nodes_reference_coords(0, 2) = 0.;
    nodes_reference_coords(1, 0) = 1.;
    nodes_reference_coords(1, 1) = 0.;
    nodes_reference_coords(1, 2) = 0.;
    nodes_reference_coords(2, 0) = 0.;
    nodes_reference_coords(2, 1) = 1.;
    nodes_reference_coords(2, 2) = 0.;
    nodes_reference_coords(3, 0) = 0.;
    nodes_reference_coords(3, 1) = 0.;
    nodes_reference_coords(3, 2) = 1.;

    return nodes_reference_coords;
}

Eigen::MatrixXf get_gauss_points()
{
    Eigen::MatrixXf gauss_points(4, 3);
    gauss_points(0, 0) = 0.5854101966;
    gauss_points(0, 1) = 0.1381966011;
    gauss_points(0, 2) = 0.1381966011;
    gauss_points(1, 0) = 0.1381966011;
    gauss_points(1, 1) = 0.5854101966;
    gauss_points(1, 2) = 0.1381966011;
    gauss_points(2, 0) = 0.1381966011;
    gauss_points(2, 1) = 0.1381966011;
    gauss_points(2, 2) = 0.5854101966;
    gauss_points(3, 0) = 0.1381966011;
    gauss_points(3, 1) = 0.1381966011;
    gauss_points(3, 2) = 0.1381966011;

    return gauss_points;
}

Eigen::VectorXf get_gauss_weights()
{
    Eigen::VectorXf gauss_weights(4);
    gauss_weights(0) = 0.0416666667;
    gauss_weights(1) = 0.0416666667;
    gauss_weights(2) = 0.0416666667;
    gauss_weights(3) = 0.0416666667;

    return gauss_weights;
}

Eigen::MatrixXf get_shapefun_coeffs()
{
    Eigen::MatrixXf nodes_reference_coords = get_reference_coords();
    float x0 = nodes_reference_coords(0, 0);
    float y0 = nodes_reference_coords(0, 1);
    float z0 = nodes_reference_coords(0, 2);
    float x1 = nodes_reference_coords(1, 0);
    float y1 = nodes_reference_coords(1, 1);
    float z1 = nodes_reference_coords(1, 2);
    float x2 = nodes_reference_coords(2, 0);
    float y2 = nodes_reference_coords(2, 1);
    float z2 = nodes_reference_coords(2, 2);
    float x3 = nodes_reference_coords(3, 0);
    float y3 = nodes_reference_coords(3, 1);
    float z3 = nodes_reference_coords(3, 2);

    Eigen::MatrixXf mat_A(4, 4);
    mat_A(0, 0) = 1.;
    mat_A(0, 1) = x0;
    mat_A(0, 2) = y0;
    mat_A(0, 3) = z0;
    mat_A(1, 0) = 1.;
    mat_A(1, 1) = x1;
    mat_A(1, 2) = y1;
    mat_A(1, 3) = z1;
    mat_A(2, 0) = 1.;
    mat_A(2, 1) = x2;
    mat_A(2, 2) = y2;
    mat_A(2, 3) = z2;
    mat_A(3, 0) = 1.;
    mat_A(3, 1) = x3;
    mat_A(3, 2) = y3;
    mat_A(3, 3) = z3;

    Eigen::MatrixXf mat_I = Eigen::MatrixXf::Identity(4, 4);

    Eigen::MatrixXf shapefun_coeffs = mat_A.ldlt().solve(mat_I);

    return shapefun_coeffs;
}

float shapefun_value(int node_idx, Eigen::Vector3f X, Eigen::MatrixXf shapefun_coeffs)
{
    float value = shapefun_coeffs(0, node_idx) 
        + shapefun_coeffs(1, node_idx) * X(0) 
        + shapefun_coeffs(2, node_idx) * X(1) 
        + shapefun_coeffs(3, node_idx) * X(2);
    
    return value;
}

// derivative_coord_idx can have value 1, 2 or 3, referring to either the x-, y- or z-derivative.
float derivative_shapefun_value(int node_idx, int derivative_coord_idx, Eigen::MatrixXf shapefun_coeffs)
{
    float value = shapefun_coeffs(derivative_coord_idx, node_idx);

    return value;
}

Eigen::MatrixXf compute_mat_Ee(Eigen::Vector3f X)
{
    Eigen::MatrixXf shapefun_coeffs = get_shapefun_coeffs();
    Eigen::Matrix3f mat_I = Eigen::Matrix3f::Identity();
    Eigen::Matrix3f mat_E0 = shapefun_value(0, X, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_E1 = shapefun_value(1, X, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_E2 = shapefun_value(2, X, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_E3 = shapefun_value(3, X, shapefun_coeffs) * mat_I;

    Eigen::MatrixXf mat_Ee(3, 12);
    std::vector<int> inds0 = {0, 1, 2};
    std::vector<int> inds1 = {3, 4, 5};
    std::vector<int> inds2 = {6, 7, 8};
    std::vector<int> inds3 = {9, 10, 11};
    mat_Ee(inds0, inds0) = mat_E0;
    mat_Ee(inds0, inds1) = mat_E1;
    mat_Ee(inds0, inds2) = mat_E2;
    mat_Ee(inds0, inds3) = mat_E3;

    return mat_Ee;
}

Eigen::MatrixXf compute_mat_De()
{
    Eigen::MatrixXf shapefun_coeffs = get_shapefun_coeffs();
    Eigen::Matrix3f mat_I = Eigen::Matrix3f::Identity();

    std::vector<int> inds0 = {0, 1, 2};
    std::vector<int> inds1 = {3, 4, 5};
    std::vector<int> inds2 = {6, 7, 8};
    std::vector<int> inds3 = {9, 10, 11};

    Eigen::Matrix3f mat_D0dx = derivative_shapefun_value(0, 1, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D1dx = derivative_shapefun_value(1, 1, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D2dx = derivative_shapefun_value(2, 1, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D3dx = derivative_shapefun_value(3, 1, shapefun_coeffs) * mat_I;

    Eigen::MatrixXf mat_Ddx(3, 12);
    mat_Ddx(inds0, inds0) = mat_D0dx;
    mat_Ddx(inds0, inds1) = mat_D1dx;
    mat_Ddx(inds0, inds2) = mat_D2dx;
    mat_Ddx(inds0, inds3) = mat_D3dx;

    Eigen::Matrix3f mat_D0dy = derivative_shapefun_value(0, 2, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D1dy = derivative_shapefun_value(1, 2, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D2dy = derivative_shapefun_value(2, 2, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D3dy = derivative_shapefun_value(3, 2, shapefun_coeffs) * mat_I;

    Eigen::MatrixXf mat_Ddy(3, 12);
    mat_Ddy(inds0, inds0) = mat_D0dy;
    mat_Ddy(inds0, inds1) = mat_D1dy;
    mat_Ddy(inds0, inds2) = mat_D2dy;
    mat_Ddy(inds0, inds3) = mat_D3dy;

    Eigen::Matrix3f mat_D0dz = derivative_shapefun_value(0, 3, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D1dz = derivative_shapefun_value(1, 3, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D2dz = derivative_shapefun_value(2, 3, shapefun_coeffs) * mat_I;
    Eigen::Matrix3f mat_D3dz = derivative_shapefun_value(3, 3, shapefun_coeffs) * mat_I;

    Eigen::MatrixXf mat_Ddz(3, 12);
    mat_Ddz(inds0, inds0) = mat_D0dz;
    mat_Ddz(inds0, inds1) = mat_D1dz;
    mat_Ddz(inds0, inds2) = mat_D2dz;
    mat_Ddz(inds0, inds3) = mat_D3dz;

    std::vector<int> inds4 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    Eigen::MatrixXf mat_D(9, 12);
    mat_D(inds0, inds4) = mat_Ddx;
    mat_D(inds1, inds4) = mat_Ddy;
    mat_D(inds2, inds4) = mat_Ddz;

    return mat_D;
}

MatsGauss compute_mats_gauss()
{
    Eigen::MatrixXf gauss_points = get_gauss_points();
    std::vector<Eigen::MatrixXf> vec_mat_Ee_gauss;
    std::vector<Eigen::MatrixXf> vec_mat_EeTEe_gauss;
    std::vector<Eigen::MatrixXf> vec_mat_De_gauss;
    std::vector<int> inds = {0, 1 ,2};

    for (size_t i = 0; i < gauss_points.rows(); i++)
    {
        Eigen::Vector3f gauss_point = gauss_points(i, inds);

        Eigen::MatrixXf mat_Ee = compute_mat_Ee(gauss_point);
        vec_mat_Ee_gauss.push_back(mat_Ee);
        Eigen::MatrixXf mat_EeTEe = mat_Ee.transpose() * mat_Ee;
        vec_mat_EeTEe_gauss.push_back(mat_EeTEe);
        Eigen::MatrixXf mat_De = compute_mat_De();
        vec_mat_De_gauss.push_back(mat_De);
    }

    MatsGauss mats_gauss = MatsGauss(vec_mat_Ee_gauss, vec_mat_EeTEe_gauss, vec_mat_De_gauss);
    
    return mats_gauss;
}

SpMat double_slice_spmat(SpMat mat, std::vector<int> row_indices, std::vector<int> col_indices)
{
    SpMat row_slice(row_indices.size(), mat.cols());
    for (size_t i = 0; i < row_indices.size(); i++)
    {
        row_slice.row(i) = mat.row(row_indices[i]);
    }
    SpMat col_slice(row_indices.size(), col_indices.size());
    for (size_t j = 0; j < col_indices.size(); j++)
    {
        col_slice.col(j) = row_slice.col(col_indices[j]);
    }
    
    return col_slice;
}