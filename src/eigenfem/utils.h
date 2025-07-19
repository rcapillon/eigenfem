//
// utils.h
//

#ifndef utils_h
#define utils_h

#include <utility>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"
#include "../../third-party/eigen-3.4.0/Eigen/Dense"


typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> RMSpMat;

Eigen::MatrixXf compute_mat_G();
Eigen::MatrixXf compute_mat_P();

Eigen::MatrixXf get_reference_coords();
Eigen::MatrixXf get_gauss_points();
Eigen::VectorXf get_gauss_weights();
Eigen::MatrixXf get_shapefun_coeffs();

float shapefun_value(int node_idx, Eigen::Vector3f X, Eigen::MatrixXf shapefun_coeffs);
float derivative_shapefun_value(int node_idx, int derivative_coord_idx, Eigen::MatrixXf shapefun_coeffs);

Eigen::MatrixXf compute_mat_Ee(Eigen::Vector3f X);
Eigen::MatrixXf compute_mat_De();

struct MatsGauss{
    std::vector<Eigen::MatrixXf> vec_mat_Ee_gauss;
    std::vector<Eigen::MatrixXf> vec_mat_EeTEe_gauss;
    std::vector<Eigen::MatrixXf> vec_mat_De_gauss;

    MatsGauss() {};
    MatsGauss(std::vector<Eigen::MatrixXf> vec_Ee, std::vector<Eigen::MatrixXf> vec_EeTEe, std::vector<Eigen::MatrixXf> vec_De)
    {
        vec_mat_Ee_gauss = vec_Ee;
        vec_mat_EeTEe_gauss = vec_EeTEe;
        vec_mat_De_gauss = vec_De;
    };
};
MatsGauss compute_mats_gauss();

SpMat double_slice_spmat(SpMat mat, std::vector<int> row_indices, std::vector<int> col_indices);

bool pair_comparator(const std::pair<float, int>& lhs, const std::pair<float, int>& rhs);

int find_index(std::vector<int>& vec, int val);

#endif