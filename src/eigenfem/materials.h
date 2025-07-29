//
// materials.h
//

// This code is subject to the terms of the MIT License.
// If a copy of the MIT License was not distributed with this file, 
// you can obtain one at https://www.mit.edu/~amini/LICENSE.md.

#ifndef materials_h
#define materials_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"


class Material
{
    public:
        Material() {};
        Material(float rho, float Y, float nu) 
        {
            Material::rho = rho;
            Material::young = Y;
            Material::poisson = nu;
            Material::mat_C = compute_mat_C();
        };
        ~Material() {};

        float rho;
        float young;
        float poisson;
        Eigen::MatrixXf mat_C;

        Eigen::MatrixXf compute_mat_C();
};

#endif