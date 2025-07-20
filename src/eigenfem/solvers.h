//
// solvers.h
//

#ifndef solvers_h
#define solvers_h

#include <iostream>
#include <utility>
#include <cmath>

#include "../../third-party/eigen-3.4.0/Eigen/Core"
#include "../../third-party/eigen-3.4.0/Eigen/SparseCholesky"
#include "../../third-party/spectra-1.1.0/include/Spectra/SymGEigsSolver.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseSymMatProd.h"
#include "../../third-party/spectra-1.1.0/include/Spectra/MatOp/SparseRegularInverse.h"

#include "model.h"
#include "utils.h"


class ModalSolver
{
    public:
        ModalSolver() {};
        ModalSolver(Model mdl);
        ~ModalSolver() {};

        Model model;

        std::vector<float> vec_freqs;
        Eigen::MatrixXf mat_modes;

        void solve(int n_modes);
};

class LinearStaticsSolver
{
    public:
        LinearStaticsSolver() {};
        LinearStaticsSolver(Model mdl);
        ~LinearStaticsSolver() {};

        Model model;

        Eigen::VectorXf U;

        void solve();
};

class FrequencySweepSolver
{
    public:
        FrequencySweepSolver() {};
        FrequencySweepSolver(Model mdl);
        ~FrequencySweepSolver() {};

        Model model;

        int n_modes;
        Eigen::MatrixXf mat_rom_basis;
        Eigen::MatrixXf mat_Mrom;
        Eigen::MatrixXf mat_Krom;
        Eigen::MatrixXf mat_Drom;
        Eigen::VectorXf vec_From;
        Eigen::MatrixXcf mat_Kdyn_rom;

        Eigen::MatrixXf mat_U_modulus;

        void load_rom(Eigen::MatrixXf rom_basis);
        void compute_rom(int n);
        void compute_Kdyn(float w);
        void solve(std::vector<float> angular_freqs);
};

#endif