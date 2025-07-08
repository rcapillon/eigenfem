# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for sparse linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

## Dependencies
- Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page
- Spectra: https://spectralib.org/

## Features
- Supported mesh types: 3D tetrahedral mesh generated using gmsh and saved in Matlab format
- Linear statics problems
- Computation of linear elastic modes and corresponding eigenfrequencies
- Frequency-domain dynamic analysis