//
// element.h
//

#ifndef element_h
#define element_h

#include <cmath>

#include "../../third-party/eigen-3.4.0/Eigen/Core"

#include "utils.h"
#include "materials.h"


struct MatsJ
{
    float det_J;
    Eigen::MatrixXf mat_invJJJ;

    MatsJ() {};
    MatsJ(float determinant_J, Eigen::MatrixXf matrix_invJJJ)
    {
        det_J = determinant_J;
        mat_invJJJ = matrix_invJJJ;
    };
};

class Element
{
    public:
        Element() {};
        Element(int num, Material mat, std::vector<int> nodes_numbers, Eigen::MatrixXf nodes_coordinates) {};
        ~Element() {};

        int number;
        Material material;
        std::vector<int> nodes_num;
        Eigen::MatrixXf nodes_coords;
        int n_nodes;
        int n_dofs;
        Eigen::MatrixXf nodes_reference_coords;
        Eigen::MatrixXf gauss_points;
        Eigen::MatrixXf gauss_weights;
        std::vector<int> dofs_num;
        Eigen::VectorXf vec_nodes_coords;
        MatsGauss mats_gauss;

        MatsJ compute_jacobian_at_gauss_points(int gauss_idx) {};
        Eigen::MatrixXf compute_mat_Me() {};
        Eigen::MatrixXf compute_mat_Ke() {};
};

#endif