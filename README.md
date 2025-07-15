# EIGENFEM
Finite Element solver written in C++ using Eigen and Spectra.
Eigen is used for linear algebra and Spectra is used to solve sparse generalized eigenvalue problems.

**IMPORTANT NOTE**: The code is still in development and not usable at this time.

## Dependencies
- Eigen: https://eigen.tuxfamily.org/
- Spectra: https://spectralib.org/

## Intended features
- Supported mesh types: 3D tetrahedral mesh generated using gmsh and saved in Matlab format
- Surface and volume forces applicable to tagged surfaces and volumes (physical groups)
- Linear statics problems
- Computation of linear elastic modes and corresponding eigenfrequencies
- Frequency-domain dynamic analysis with a reduced-order model using elastic modes

## Limitations
- Only fully tetrahedral meshes are handled
- Nodal forces are not handled because the code doesn't look for 0D physical groups in the mesh file
- Surface and volume forces must be constant vectors for each 2D or 3D physical groups
- Only one 3D physical group is handled, which forces to use a single material for the whole domain

## Documentation and usage
Documentation and user guide will be put in the wiki of the github repository.