//
// element.cpp
//

#include "element.h"


Eigen::MatrixXf mat_G = compute_mat_G();
Eigen::MatrixXf mat_P = compute_mat_P();
Eigen::MatrixXf gauss_points = get_gauss_points();
Eigen::MatrixXf gauss_weights = get_gauss_weights();
Eigen::MatrixXf nodes_reference_coords = get_reference_coords();
Eigen::MatrixXf shapefun_coeffs = get_shapefun_coeffs();
