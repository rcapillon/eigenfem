# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

**IMPORTANT NOTE**: 
The code is still in development and not fully usable at this time. Only specific demonstration cases can be run. 
Installation procedure will be specified once all 4 example cases can be compiled and run without bugs.
Usage documentation will be added later when custom computations can be performed using a custom input file format.

## Dependencies
- Eigen v3.4.0: https://eigen.tuxfamily.org/
    - Licensed under MPL2
- Spectra v1.1.0: https://spectralib.org/
    - Licensed under MPL2
    - Copyright 2015-2025, Yixuan Qiu

## Intended features
- Supported mesh type: 3D tetrahedral mesh generated using gmsh and saved in Matlab format
- Surface and volume forces applicable to tagged surfaces and volumes (physical groups)
- Linear statics problems
- Computation of linear elastic modes and corresponding eigenfrequencies
- Frequency-domain dynamic analysis with a reduced-order model using elastic modes
- Output deformed meshes as VTK files (recommended software for visualization is Paraview: https://www.paraview.org/)

## Limitations
- Only fully tetrahedral (first order) meshes are handled
- No multi-threading, parallelization of any sort
- Only zero Dirichlet boundary conditions are handled
- Nodal forces are not handled because the code doesn't look for 0D physical groups in the mesh file
- Only one 3D physical group is handled, which forces to use a single material for the whole domain

## Documentation and usage
Documentation and user guide will be put in the wiki of the github repository.