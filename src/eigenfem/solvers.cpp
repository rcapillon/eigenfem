//
// solvers.cpp
//

#include <iostream>
#include <utility>
#include <cmath>

#include "solvers.h"

// #include "../../third-party/eigen-3.4.0/Eigen/Eigenvalues"


ModalSolver::ModalSolver(Model mdl)
{
    model = mdl;
}

void ModalSolver::solve_test()
{
    model.create_dof_lists();
    model.assemble_M_K();
    model.apply_dirichlet();

    Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> ges;
    ges.compute(Eigen::MatrixXf(model.mat_Kff), Eigen::MatrixXf(model.mat_Mff));
    Eigen::VectorXcf vec_eigvals = ges.eigenvalues();
    test_freqs = Eigen::VectorXcf(vec_eigvals.size());
    for (size_t i = 0; i < vec_eigvals.size(); i++)
    {
        // test_freqs(i) = sqrt(vec_eigvals(i));
        test_freqs(i) = vec_eigvals(i);
    }
}

void ModalSolver::solve(int n_modes)
{
    model.create_dof_lists();
    model.assemble_M_K();
    model.apply_dirichlet();

    Spectra::SparseSymMatProd<float> Aop(model.mat_Kff);
    Spectra::SparseRegularInverse<float> Bop(model.mat_Mff);

    int ncv = 2 * n_modes;
    Spectra::SymGEigsSolver<
        Spectra::SparseSymMatProd<float>, 
        Spectra::SparseRegularInverse<float>, 
        Spectra::GEigsMode::RegularInverse> geigs(Aop, Bop, n_modes, ncv);
    
    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::SmallestMagn);

    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        Eigen::VectorXf eigenvalues = geigs.eigenvalues();
        mat_modes = geigs.eigenvectors();

        std::vector<float> vec_eigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
        std::vector<std::pair<float, int>> vec_eigvals_to_sort;
        for (int i = 0; i < vec_eigenvalues.size(); i++)
        {
            vec_eigvals_to_sort.push_back(std::make_pair(vec_eigenvalues[i], i));
        }
        std::sort(vec_eigvals_to_sort.begin(), vec_eigvals_to_sort.end(), pair_comparator);
        std::vector<int> sort_indices;
        for (size_t i = 0; i < vec_eigvals_to_sort.size(); i++)
        {
            vec_freqs.push_back(sqrt(vec_eigvals_to_sort[i].first));
            sort_indices.push_back(vec_eigvals_to_sort[i].second);
        }
        mat_modes = mat_modes(Eigen::all, sort_indices);
    }
    else
    {
        std::cout << "Eigensolver did not converge." << std::endl;
    }
}