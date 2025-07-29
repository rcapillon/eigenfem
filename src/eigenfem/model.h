//
// model.h
//

// This code is subject to the terms of the MIT License.
// If a copy of the MIT License was not distributed with this file, 
// you can obtain one at https://www.mit.edu/~amini/LICENSE.md.

#ifndef model_h
#define model_h

#include <tuple>
#include <cmath>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"

#include "mesh.h"
#include "utils.h"


typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::Triplet<float> triplet;

class Model
{
    public:
        Model() {};
        Model(Mesh msh, std::vector<int> dir_tags);
        Model(Mesh msh, std::vector<int> dir_tags, float a_M, float a_K);
        Model(Mesh msh, 
            std::vector<int> dir_tags, 
            std::vector<std::tuple<int, Eigen::VectorXf>> nod_forces,
            std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces,
            std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces,
            float a_M = 0., float a_K = 0.);
        ~Model() {};

        Mesh mesh;
        std::vector<int> dirichlet_tags;
        std::vector<std::tuple<int, Eigen::VectorXf>> tuples_nodal_forces;
        std::vector<std::tuple<int, Eigen::VectorXf>> tuples_surface_forces;
        std::vector<std::tuple<int, Eigen::VectorXf>> tuples_volume_forces;
        float alpha_M;
        float alpha_K;

        SpMat mat_M;
        SpMat mat_K;
        SpMat mat_D;
        Eigen::VectorXf vec_Fn;
        Eigen::VectorXf vec_Fs;
        Eigen::VectorXf vec_Fv;
        Eigen::VectorXf vec_F;

        std::vector<int> free_dofs;
        std::vector<int> dirichlet_dofs;

        SpMat mat_Mff;
        SpMat mat_Kff;
        SpMat mat_Dff;
        Eigen::VectorXf vec_Ff;

        void create_dof_lists();
        void assemble_M();
        void assemble_K();
        void assemble_M_K();
        void compute_D_Rayleigh();

        void assemble_Fn();
        void assemble_Fs();
        void assemble_Fv();
        void compute_F();

        void apply_dirichlet();
};

#endif