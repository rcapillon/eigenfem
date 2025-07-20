//
// solvers.cpp
//

#include "solvers.h"


const float PI = 3.141592653589793;

//
// MODAL SOLVER
//

ModalSolver::ModalSolver(Model mdl)
{
    model = mdl;
}

void ModalSolver::solve(int n_modes)
{
    model.create_dof_lists();
    model.assemble_M_K();
    model.apply_dirichlet();

    Spectra::SparseSymMatProd<float> Aop(model.mat_Kff);
    Spectra::SparseRegularInverse<float> Bop(model.mat_Mff);

    // very experimental scheme for choosing ncv
    int ncv = std::max(20, 4 * n_modes);
    while (ncv > model.mat_Mff.rows())
    {
        ncv /= 2;
        while (ncv > n_modes)
        {
            ncv += 1;
        }
    }
    
    Spectra::SymGEigsSolver<
        Spectra::SparseSymMatProd<float>, 
        Spectra::SparseRegularInverse<float>, 
        Spectra::GEigsMode::RegularInverse> geigs(Aop, Bop, n_modes, ncv);
    
    geigs.init();
    int nconv = geigs.compute(Spectra::SortRule::SmallestMagn);

    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        Eigen::VectorXf eigenvalues = geigs.eigenvalues();
        Eigen::MatrixXf mat_modes_free = geigs.eigenvectors();

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
            vec_freqs.push_back(sqrt(vec_eigvals_to_sort[i].first) / (2 * PI));
            sort_indices.push_back(vec_eigvals_to_sort[i].second);
        }
        mat_modes_free = mat_modes_free(Eigen::all, sort_indices);

        mat_modes = Eigen::MatrixXf::Zero(model.mesh.n_dofs, mat_modes_free.cols());
        mat_modes(model.free_dofs, Eigen::all) = mat_modes_free;
    }
    else
    {
        std::cout << "Eigensolver failed." << std::endl;
    }
}

//
// LINEAR STATICS SOLVER
//

LinearStaticsSolver::LinearStaticsSolver(Model mdl)
{
    model = mdl;
}

void LinearStaticsSolver::solve()
{
    model.create_dof_lists();
    model.assemble_K();
    model.assemble_Fs();
    model.assemble_Fv();
    model.compute_F();
    model.apply_dirichlet();

    Eigen::SimplicialLDLT<SpMat> solver;
    solver.compute(model.mat_Kff);
    if (solver.info() == Eigen::Success)
    {
        Eigen::VectorXf U_free = solver.solve(model.vec_Ff);
        if (solver.info() == Eigen::Success)
        {
            U = Eigen::VectorXf::Zero(model.mesh.n_dofs);
            U(model.free_dofs) = U_free;
        }
        else
        {
            std::cout << "Solver failed." << std::endl;
        }
    }
    else
    {
        std::cout << "Matrix decomposition failed." << std::endl;
    }
}

//
// FREQUENCY SWEEP SOLVER
//

FrequencySweepSolver::FrequencySweepSolver(Model mdl)
{
    model = mdl;
}

void FrequencySweepSolver::load_rom(Eigen::MatrixXf rom_basis)
{
    mat_rom_basis = rom_basis;
    n_modes = mat_rom_basis.cols();

    model.create_dof_lists();
    model.assemble_M_K();
    model.compute_D_Rayleigh();
    model.assemble_Fs();
    model.assemble_Fv();
    model.compute_F();
    model.apply_dirichlet();

    Eigen::MatrixXf mat_Mrom = mat_rom_basis.transpose() * model.mat_Mff * mat_rom_basis;
    Eigen::MatrixXf mat_Krom = mat_rom_basis.transpose() * model.mat_Kff * mat_rom_basis;
    Eigen::MatrixXf mat_Drom = mat_rom_basis.transpose() * model.mat_Dff * mat_rom_basis;
    Eigen::VectorXf vec_From = mat_rom_basis.transpose() * model.vec_Ff;
}

void FrequencySweepSolver::compute_rom(int n)
{
    n_modes = n;
    ModalSolver modal_solver(model);
    modal_solver.solve(n_modes);
    mat_rom_basis = modal_solver.mat_modes;

    model.create_dof_lists();
    model.assemble_M_K();
    model.compute_D_Rayleigh();
    model.assemble_Fs();
    model.assemble_Fv();
    model.compute_F();
    model.apply_dirichlet();

    Eigen::MatrixXf mat_Mrom = mat_rom_basis.transpose() * model.mat_Mff * mat_rom_basis;
    Eigen::MatrixXf mat_Krom = mat_rom_basis.transpose() * model.mat_Kff * mat_rom_basis;
    Eigen::MatrixXf mat_Drom = mat_rom_basis.transpose() * model.mat_Dff * mat_rom_basis;
    Eigen::VectorXf vec_From = mat_rom_basis.transpose() * model.vec_Ff;
}

void FrequencySweepSolver::compute_Kdyn(float w)
{
    std::complex<float> pure_imag(0, 1);
    mat_Kdyn_rom = -pow(w, 2.) * mat_Mrom + pure_imag * mat_Drom + mat_Krom;
}

void FrequencySweepSolver::solve(std::vector<float> angular_freqs)
{
    mat_U_modulus = Eigen::MatrixXf::Zero(model.mesh.n_dofs, angular_freqs.size());

    for (size_t i = 0; i < angular_freqs.size(); i++)
    {
        compute_Kdyn(angular_freqs[i]);

        Eigen::SimplicialLDLT<Eigen::MatrixXcf> solver;
        solver.compute(mat_Kdyn_rom);
        if (solver.info() == Eigen::Success)
        {
            Eigen::VectorXcf U_free = solver.solve(vec_From);
            if (solver.info() == Eigen::Success)
            {
                mat_U_modulus(model.free_dofs, i) = U_free.cwiseAbs();
            }
            else
            {
                std::cout << "Solver failed." << std::endl;
            }
        }
        else
        {
            std::cout << "Matrix decomposition failed." << std::endl;
        }
    }
    
}