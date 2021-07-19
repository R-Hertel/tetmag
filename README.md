![Github All Releases](https://img.shields.io/github/downloads/R-Hertel/tetmag/total?style=plastic)
![GitHub Forks](https://img.shields.io/github/forks/R-Hertel/tetmag?style=plastic)
![GitHub Stars](https://img.shields.io/github/stars/R-Hertel/tetmag?style=plastic)
![GitHub Release Date](https://img.shields.io/github/release-date/R-Hertel/tetmag?style=plastic)


 <img src="https://github.com/R-Hertel/tetmag/blob/main/resources/tetmagLogo_v5.png" width="400" >

`tetmag` is a finite-element software for large-scale micromagnetic simulations.
<!--- ![tetmag logo](https://github.com/R-Hertel/tetmag/blob/main/resources/tetmagLogo_v1.png) --->

       

## Installation

### Requirements

- Ubuntu 20.04: 
  Install the following


````bash 
sudo apt-get install libboost-all-dev libeigen3-dev libnetcdf-dev build-essential cmake 
```` 

- Optional:
  - GPU-acceleration with [CUDA](https://developer.nvidia.com/cuda-downloads)  (version 10.1 or higher needed)
  - Multi-threading with [OpenMP](https://www.openmp.org/):
    ```` sudo apt-get install libomp-dev ````
 



### Compilation

````bash 
 wget https://github.com/R-Hertel/tetmag/tetmag.git 
 cd tetmag && mkdir build && cd build 
 cmake ..
 make -j$(nproc)
````

- Notes:
    - Installation on other platforms and on other linux distributions should be possible, but hasn't been tested. 
    - An internet connection is required during the build process.

## Additional software
Although `tetmag` is a standalone software that -unlike plugin-type applications- does not require a specific framework to operate, its workflow involves a series of processes that must be handled by other applications. This concerns two categories of operations:

- **Preprocessing** - Definition of the sample geometry and finite-element mesh generation
- **Postprocessing** - Visualization and analysis of the results 

The following software can accomplish these tasks:

- [ParaView](https://www.paraview.org) - Data visualization 
- [gmsh](https://gmsh.info) or [netgen](https://ngsolve.org) - Finite-element mesh generation
- [FreeCAD](https://www.freecadweb.org) - Design of three-dimensional objects

The documentation describes the required pre- and postprocessing steps and shows how `tetmag` interacts with this external software. 

## Usage
 - A user guide is being prepared, and a link will be posted here once it is available. 
 - See the simulation studies in the [examples directory](https://github.com/R-Hertel/tetmag/tree/main/resources) for a quick introduction to the usage of tetmag.


## Credits
The tetmag software is powered by
- [AMGCL](https://github.com/ddemidov/amgcl) - Iterative solution of sparse systems
- [SUNDIALS / CVODE](https://github.com/LLNL/sundials) - Solution of ordinary differential equation systems
- [Eigen](https://gitlab.com/libeigen/eigen) - Linear algebra
- [H2Lib](https://github.com/H2Lib) - Hierarchical matrix data compression
- [VTK](https://github.com/Kitware/VTK) and [gmsh](https://gitlab.onelab.info/gmsh/gmsh) - I/O of finite-element mesh data

- [CUDA / Thrust](https://developer.nvidia.com/cuda-downloads) - GPU acceleration
- [boost](https://github.com/boostorg) - Powerful C++ libraries 
- [CMake](https://github.com/Kitware/CMake) - Build system generator

