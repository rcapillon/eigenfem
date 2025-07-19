# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

**IMPORTANT NOTE**: The code is still in development and not usable at this time.

## Dependencies
- Eigen: https://eigen.tuxfamily.org/
- Spectra: https://spectralib.org/

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
- Surface and volume forces must be constant vectors for each 2D or 3D physical groups
- Only one 3D physical group is handled, which forces to use a single material for the whole domain
- While the mesh is limited to a single element type, much performance is lost by computing element matrices every time like if all element types were unknown and possibly producing more variations that the 4-node tetrahedron

## Documentation and usage
Documentation and user guide will be put in the wiki of the github repository.