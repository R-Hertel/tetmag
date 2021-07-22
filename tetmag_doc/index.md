# Welcome to the tetmag documentation

`tetmag` is a standalone finite-element software for three-dimensional micromagnetic simulations. The software package allows one to simulate -within the framework of micromagnetic theory- the structure and the dynamics of the magnetization in ferromagnetic nano- and microstructures of arbitrary shape.

On these length scales, the shape of the samples can have a strong impact on their magnetic properties. The flexibility of the finite-element method allows one to accurately consider details of the geometry  -such as, e.g., surface curvature- and thereby study how these features affect the static and dynamic properties of the magnetization. Owing to a combination of the finite-element method and the boundary-element method, `tetmag` is also particularly efficient in simulating arrays of interacting magnetic particles.

The software is designed to solve large-scale problems by exploiting (optional) GPU acceleration, OpenMP parallelization, and an efficient data compression with hierarchical matrices. 

The source code is available under GNU-GPL3 License from the [tetmag repository](https://github.com/R-Hertel/tetmag) on GitHub.

## Referencing



|  |  | 
|---|---|
| [[Her19](bib/hertel09.bib)] | *Large-scale magnetostatic field calculation in finite element micromagnetics with H2-matrices*<br>R. Hertel, S. Christophersen, and S. Boerm <br>[Journal of Magnetism and Magnetic Materials 477, 118 (2019)](https://doi.org/10.1016/j.jmmm.2018.12.103).|
| [[Her21](bib/hertel16.bib)] |  *tetmag - a general-purpose micromagnetic finite-element software for large-scale simulations*<br>Riccardo Hertel<br>(2021, to be published) |
| [[Her21](bib/hertel16.bib)] |  *Micromagnetic modeling with tetmag:  an  open-source finite-element software with GPU acceleration*<br>Riccardo Hertel<br>(2021, to be published) |



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
