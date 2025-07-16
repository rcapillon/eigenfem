//
// model.h
//

#ifndef model_h
#define model_h

#include <tuple>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/Sparse"

#include "mesh.h"


class Model
{
    public:
        Model() {};
        Model(Mesh msh, 
            std::vector<int> dir_tags, 
            std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces,
            std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces,
            float a_M, float a_K) {};
        ~Model() {};

        Mesh mesh;
        std::vector<int> dirichlet_tags;
        std::vector<std::tuple<int, Eigen::VectorXf>> surface_forces;
        std::vector<std::tuple<int, Eigen::VectorXf>> volume_forces;
        float alpha_M;
        float alpha_K;
};

#endif