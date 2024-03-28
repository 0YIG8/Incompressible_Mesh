# Incompressible_Mesh

This repository contains all the code files for the SIMPLE Algorithm Assignment.

The code generates a computational mesh using input data about points, faces, cells, and boundaries generated using OpenFOAM. It is capable of discretizing terms like convection, diffusion (Laplacian), and temporal derivatives, thus being able to solve the incompressible Navier-Stokes equations. The Eigen C++ library has been used to store and solve the matrices used in the ScalarField and VectorField classes.

The test_main.C file contains the SIMPLE algorithm for pressure-velocity coupling for a simple 40x40 2D lid driven cavity test case. Input files for other mesh resolutions are also available in the other folders provided. The code is compiled using CMake; but can be run alternatively on the Windows terminal using the bash script provided. For the bash script, by default, the code runs with debug mode ON -> need to comment out the debug mode commands and uncomment the others.

For the CMake option, to compile and build, run:
```bash
$ cmake --build build
$ ./build/main
```

For the bash script, to compile and build, run:
```cmd
./run
```
