# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

**IMPORTANT NOTE**: 
The code has not been validated on test cases. Demonstration cases can be run and custom simulations can be performed, but usage documentation will be added soon to the repository's wiki.

## Dependencies
- Eigen v3.4.0: https://eigen.tuxfamily.org/
    - Licensed under MPL2
- Spectra v1.1.0: https://spectralib.org/
    - Licensed under MPL2
    - Copyright 2015-2025, Yixuan Qiu

## Installation
Requirements: 
- CMake (minimal version required: 3.28.3)
- OS: MacOS (tested on MacOS Sequoia 15.5 with a M1 2020 Macbook Pro), Linux (tested on Ubuntu arm64 24.04.2).

No testing on Windows has been done for now.

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

## Running example scripts
The installation produces 6 executables in the build folder:
- eigenfem: requires custom input files to run a simulation
- example_modal, example_statics, example_statics_nodal, example_frequency_coarse, example_frequency_fine: examples showcasing the code's features and the 3 available solvers.

You can run one of the examples from the build folder with, for instance:
```
./example_modal
```
A short description of these examples can be found in the repository's Wiki.
Results (VTK files of deformed meshes) will be created in subfolders of the /examples/results/ folder.
The recommended software to visualize output VTK files is Paraview: https://www.paraview.org/

## Running custom simulations
Running a custom simulation requires using the eigenfem executable in the build folder and passing it a text input file as argument. 5 example input files are given in the /examples/ folder (.input files), performing simulations very close to the example scripts described in the previous section.

Documentation on how to write a custom input file can be found on the repository's [wiki](https://github.com/rcapillon/eigenfem/wiki/Custom-simulation). Users can also read the sample input files to understand how to create their own.

In order to run a custom simulation, assuming the eigenfem executable is in the /build/ folder and the input file is located in the examples folder and named 'custom.input' (the file's extension is irrelevant), run the following command inside the /build/ folder:
```
./eigenfem ../examples/custom.input
```

To run a simulation with one of the pre-made input files, go inside the /build/ folder and for instance use:
```
./eigenfem ../examples/example_input_modal.input
```

The input file specifies a directory to store VTK output files. In the provided examples, this folder is set as '/examples/results/example_eigenfem/'.

## Current features
- Supported mesh type: 3D tetrahedral mesh generated using gmsh and saved in Matlab format
- Nodal, surface and volume forces applicable to tagged surfaces and volumes (physical groups)
- Linear statics problems
- Computation of linear elastic modes and corresponding eigenfrequencies
- Frequency-domain dynamic analysis with a reduced-order model using elastic modes (or a loaded reduced-order basis)
- Saving or loading a basis of modes to use for reduced-order modeling with the frequency-domain solver
- Output deformed meshes as VTK files

## Limitations
- Only fully tetrahedral (first order) meshes are handled
- No multi-threading, parallelization of any sort
- Only zero Dirichlet boundary conditions are handled
- Only one 3D physical group is handled, which forces to use a single material for the whole domain
- Only isotropic elastic materials are handled

## (Non-exhaustive) list of intended future features
- Handling of multiple 3D physical groups to enable domains with multiple materials or to isolate 3D elements to, for instance, calculate stress
- Handling of first-order hexahedral (6 rectangular faces) and prism (2 triangular faces and 3 rectangular faces) elements