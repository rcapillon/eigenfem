//
// solvers.h
//

#ifndef solvers_h
#define solvers_h

#include "../../third-party/eigen-3.4.0/Eigen/Core"
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

#endif