//
// materials.h
//

#ifndef materials_h
#define materials_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"


class Material
{
    public:
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