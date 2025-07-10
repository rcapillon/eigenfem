//
// utils.h
//

#ifndef utils_h
#define utils_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/Dense"


Eigen::MatrixXf compute_mat_G();
Eigen::MatrixXf compute_mat_P();
Eigen::MatrixXf get_reference_coords();
Eigen::MatrixXf get_gauss_points();
Eigen::MatrixXf get_gauss_weights();
Eigen::MatrixXf get_shapefun_coeffs();
float shapefun_value(int node_idx, Eigen::Vector3f X, Eigen::MatrixXf shapefun_coeffs);
float derivative_shapefun_value(int node_idx, int derivative_coord_idx, Eigen::MatrixXf shapefun_coeffs);
Eigen::MatrixXf compute_mat_Ee(Eigen::Vector3f X);

#endif