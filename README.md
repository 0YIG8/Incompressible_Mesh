# Incompressible_Mesh

This repository contains all the code files for the SIMPLE Algorithm Assignment.

The code generates a computational mesh using input data about points, faces, cells, and boundaries, and should work for unstructured meshes. It is capable of discretizing terms like convection, diffusion (Laplacian), and temporal derivatives, thus being able to solve the incompressible Navier-Stokes equations. 

The test_main.C file contains the SIMPLE algorithm for pressure-velocity coupling for a simple 40x40 lid driven cavity test case.

To compile, run:
rm -rf build
cmake -S . -B build
cmake --build build
./build/main