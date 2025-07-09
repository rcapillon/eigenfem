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

Eigen::MatrixXf get_gauss_weights()
{
    Eigen::VectorXf gauss_weights(4);
    gauss_weights(0) = 0.0416666667;
    gauss_weights(1) = 0.0416666667;
    gauss_weights(2) = 0.0416666667;
    gauss_weights(3) = 0.0416666667;

    return gauss_weights;
}
