![Github All Releases](https://img.shields.io/github/downloads/R-Hertel/tetmag/total?style=plastic)
![GitHub Forks](https://img.shields.io/github/forks/R-Hertel/tetmag?style=plastic)
![GitHub Stars](https://img.shields.io/github/stars/R-Hertel/tetmag?style=plastic)
![GitHub Release Date](https://img.shields.io/github/release-date/R-Hertel/tetmag?style=plastic)

# tetmag
A finite-element software for large-scale micromagnetic simulations.
<!--- ![tetmag logo](https://github.com/R-Hertel/tetmag/blob/main/resources/tetmagLogo_v1.png) --->
 <img src="https://github.com/R-Hertel/tetmag/blob/main/resources/tetmagLogo_v1.png" width="400"  >
       

## Install

### Requirements:

- Ubuntu 20.04: 
  Install the following

````bash 
sudo apt-get install libboost-all-dev libeigen3-dev libnetcdf-dev build-essential cmake 
```` 

- Optional:
  - GPU-acceleration with [CUDA](https://developer.nvidia.com/cuda-downloads)  > 10.1
  - Multi-threading with OpenMP 
    ```` sudo apt-get install libomp-dev ````
 



### Compilation:

````bash 
 wget https://github.com/R-Hertel/tetmag/tetmag.git 
 cd tetmag && mkdir build && cd build 
 cmake ..
 make -j$(nproc)
````
## Pre- and Postprocessing:

- [ParaView](https://www.paraview.org) - Data visualization 
- [gmsh](https://www.paraview.org/) or [netgen](https://ngsolve.org) - Finite-element mesh generation
- [FreeCAD](https://www.freecadweb.org) - Design of three-dimensional objects

## Notes
Compilation on other platforms and on other linux distributions should be possible, but hasn't been tested. 


## Credits
The tetmag software is powered by:
- [AMGCL](https://github.com/ddemidov/amgcl) - Iterative solution of sparse systems
- [SUNDIALS / CVODE](https://github.com/LLNL/sundials) - Solution of ordinary differential equation (ODE) systems
- [Eigen](https://gitlab.com/libeigen/eigen) - Linear algebra
- [H2Lib](https://github.com/H2Lib) - Hierarchical matrix data compression
- [VTK](https://github.com/Kitware/VTK) and [gmsh](https://gitlab.onelab.info/gmsh/gmsh) - I/O of finite-element mesh data

- [CUDA / Thrust](https://developer.nvidia.com/cuda-downloads) - GPU acceleration
- [boost](https://github.com/boostorg) - Powerful C++ libraries 
- [CMake](https://github.com/Kitware/CMake) - Build system generator
