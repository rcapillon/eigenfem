# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

**IMPORTANT NOTE**: 
The code is still in development and not fully usable at this time. The code has not been validated on test cases. For now, only specific demonstration cases can be run. 
Custom usage documentation will be added later.

## Dependencies
- Eigen v3.4.0: https://eigen.tuxfamily.org/
    - Licensed under MPL2
- Spectra v1.1.0: https://spectralib.org/
    - Licensed under MPL2
    - Copyright 2015-2025, Yixuan Qiu

## Installation and running example scripts
Installation can be performed using CMake (minimal version required: 3.28.3). 
The installation procedure has been tested on MacOS Sequoia 15.5 and Ubuntu 24.04.2, but should work for other versions as well. No testing on Windows has been done for now.


First, clone the repository in the folder of your choice with:
```
git clone git@github.com:rcapillon/eigenfem.git
```
Then, use the following commands:
```
cd eigenfem/
mkdir build
cd build/
cmake ..
cmake --build .
```
This will produce 5 executables in the build folder:
- eigenfem: meant to use custom input files to run a custom simulation using one of the 3 available solvers (will be available later)
- example_modal, example_statics, example_frequency_coarse, example_frequency_fine: examples showcasing the code's features and the 3 available solvers.

You can then run one of the examples from the build folder with, for instance:
```
./example_modal
```
A short description of these examples can be found in the repository's Wiki.
Results (VTK files of deformed meshes) will be created in subfolders of the /examples/results folder.
The recommended software to visualize output VTK files is Paraview: https://www.paraview.org/

## Current features
- Supported mesh type: 3D tetrahedral mesh generated using gmsh and saved in Matlab format
- Surface and volume forces applicable to tagged surfaces and volumes (physical groups)
- Linear statics problems
- Computation of linear elastic modes and corresponding eigenfrequencies
- Frequency-domain dynamic analysis with a reduced-order model using elastic modes (or a loaded reduced-order basis)
- Saving or loading a basis of modes to use for reduced-order modeling with the frequency-domain solver
- Output deformed meshes as VTK files

## Limitations
- Only fully tetrahedral (first order) meshes are handled
- No multi-threading, parallelization of any sort
- Only zero Dirichlet boundary conditions are handled
- Nodal forces are not handled because the code doesn't look for 0D physical groups in the mesh file
- Only one 3D physical group is handled, which forces to use a single material for the whole domain
- Only isotropic elastic materials are handled
- It is not currently possible to define and run a fully custom simulation using one of the available solvers

## (Non-exhaustive) list of intended future features
- Handling of 0D physical groups to enable nodal forces and plotting of nodal displacement values
- Handling of multiple 3D physical groups to enable domains with multiple materials or to isolate 3D elements to, for instance, calculate stress
- Handling of first-order hexahedral (6 rectangular faces) and prism (2 triangular faces and 3 rectangular faces) elements