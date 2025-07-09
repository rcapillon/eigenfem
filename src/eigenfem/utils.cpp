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
