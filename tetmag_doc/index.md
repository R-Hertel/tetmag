# Welcome to the tetmag documentation

`tetmag` is a standalone finite-element software for three-dimensional micromagnetic simulations. The software package allows one to simulate -within the framework of micromagnetic theory- the structure and the dynamics of the magnetization in ferromagnetic nano- and microstructures of arbitrary shape.

On these length scales, the shape of the samples can have a strong impact on their magnetic properties. The flexibility of the finite-element method allows one to accurately consider details of the geometry  -such as, e.g., surface curvature- and thereby study how these features affect the static and dynamic properties of the magnetization. Owing to a combination of the finite-element method and the boundary-element method, `tetmag` is also efficient in simulating arrays of interacting magnetic particles.

The software is designed to solve large-scale problems by exploiting (optional) GPU acceleration, OpenMP parallelization, and an efficient data compression with hierarchical matrices. 

The source code is available under AGPLv3 License from the [tetmag repository](https://github.com/R-Hertel/tetmag) on GitHub.


## Credits


`tetmag` is powered by these awesome software projects:

- [Eigen](https://gitlab.com/libeigen/eigen) - Linear algebra
- [H2Lib](https://github.com/H2Lib) - Hierarchical matrix data compression
- [AMGCL](https://github.com/ddemidov/amgcl) - Iterative solution of sparse systems
- [SUNDIALS / CVODE](https://github.com/LLNL/sundials) - Solution of ordinary differential equation systems
- [VTK](https://github.com/Kitware/VTK) and [gmsh](https://gitlab.onelab.info/gmsh/gmsh) - I/O routines for finite-element mesh data

- [CUDA / Thrust](https://developer.nvidia.com/cuda-downloads) - GPU acceleration
- [boost](https://github.com/boostorg) - Powerful C++ libraries 
- [CMake](https://github.com/Kitware/CMake) - Build system generator
