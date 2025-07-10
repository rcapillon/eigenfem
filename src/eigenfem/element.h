//
// element.h
//

#ifndef element_h
#define element_h

#include <cmath>

#include "../../third-party/eigen-3.4.0/Eigen/Core"

#include "utils.h"
#include "materials.h"


Eigen::MatrixXf mat_G;
Eigen::MatrixXf mat_P;
Eigen::MatrixXf gauss_points;
Eigen::MatrixXf gauss_weights;
Eigen::MatrixXf nodes_reference_coords;
Eigen::MatrixXf shapefun_coeffs;

MatsGauss mats_gauss;

class Element
{
    public:
        Element(int number, Material material, std::vector<int> nodes_num, Eigen::Vector3f nodes_coords)
        {
        };
        ~Element() {};

};

#endif