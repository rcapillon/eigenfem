//
// utils.h
//

#ifndef utils_h
#define utils_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"


Eigen::MatrixXf compute_mat_G();
Eigen::MatrixXf compute_mat_P();
Eigen::MatrixXf get_reference_coords();
Eigen::MatrixXf get_gauss_points();
Eigen::MatrixXf get_gauss_weights();

#endif