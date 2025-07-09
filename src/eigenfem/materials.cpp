//
// materials.cpp
//

#include "materials.h"


 Eigen::MatrixXf Material::compute_mat_C()
{
    float Y = Material::young;
    float nu = Material::poisson;
    float lame1 = (Y * nu) / ((1 + nu) * (1 - 2 * nu));
    float lame2 = Y / (2 * (1 + nu));

    Eigen::MatrixXf mat_C(6, 6);
    mat_C(0, 0) = lame1;
    mat_C(0, 1) = lame1;
    mat_C(0, 2) = lame1;
    mat_C(1, 0) = lame1;
    mat_C(1, 1) = lame1;
    mat_C(1, 2) = lame1;
    mat_C(2, 0) = lame1;
    mat_C(2, 1) = lame1;
    mat_C(2, 2) = lame1;
    mat_C(0, 0) += 2 * lame2;
    mat_C(1, 1) += 2 * lame2;
    mat_C(2, 2) += 2 * lame2;
    mat_C(3, 3) = 2 * lame2;
    mat_C(4, 4) = 2 * lame2;
    mat_C(5, 5) = 2 * lame2;

    return mat_C;
}
