//
// model.h
//

#ifndef model_h
#define model_h

#include <tuple>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCore"

#include "mesh.h"


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
            std::vector<std::tuple<int, Eigen::VectorXf>> surf_forces,
            std::vector<std::tuple<int, Eigen::VectorXf>> vol_forces,
            float a_M, float a_K);
        ~Model() {};

        Mesh mesh;
        std::vector<int> dirichlet_tags;
        std::vector<std::tuple<int, Eigen::VectorXf>> surface_forces;
        std::vector<std::tuple<int, Eigen::VectorXf>> volume_forces;
        float alpha_M;
        float alpha_K;

        SpMat mat_M;
        SpMat mat_K;
        SpMat mat_D;

        std::vector<int> free_dofs;
        std::vector<int> dirichlet_dofs;

        SpMat mat_Mff;
        SpMat mat_Kff;
        SpMat mat_Dff;

        void create_dof_lists();
        void assemble_M();
        void assemble_K();
        void assemble_M_K();
        void compute_D_Rayleigh();
        void apply_dirichlet();
};

#endif