# Installation notes for tetmag

This page describes how to compile `tetmag`. The installation guidelines described here refer to **ubuntu 20.04 LTS**. It should be possible to compile `tetmag` also on other Linux distributions and on Windows (64bit) systems. However, `tetmag` has not yet been tested on systems other than ubuntu and MacOS.  

## Prerequisites

As a minimal requirement, the following libraries and utilities should be installed:

    sudo apt-get install libboost-all-dev libeigen3-dev libnetcdf-dev libvtk7-dev build-essential cmake 

## Standard installation 

Once the software listed above is installed, `tetmag` can be compiled with

    git clone https://github.com/R-Hertel/tetmag.git
    cd tetmag && mkdir build && cd build 
    cmake ..
    make -j$(nproc)

The complation should end after some time with the message

````bash
...
[100%] Linking CXX executable tetmag
[100%] Built target tetmag
````

!!! note
    Several external libraries will be downloaded and installed during the compilation. An active internet connection is therefore required to compile the code.

The executable generated with this standard procedure can be used for CPU-based micromagnetic simulations with OpenMP parallelization. 
A few simple cases showing how to run simulations with `tetmag` are discussed in the [examples section](examples.md).

## Installation with CUDA 

The `tetmag` software can exploit GPU acceleration based on NVIDIA's [CUDA toolkit](https://developer.nvidia.com/cuda-toolkit). This feature requires a CUDA-compatible graphics card and an installation of CUDA (version 10.2 or higher). Please follow the instructions on the [NVIDIA developer page](https://developer.nvidia.com/cuda-downloads) on how to install the CUDA Toolkit and reboot the computer after completing the CUDA installation.

To compile a CUDA-accelerated version of `tetmag`, the procedure is as described above, except for an additional compiler flag `-DUSE_CUDA` passed to `cmake`:

    cd tetmag && mkdir build && cd build 
    cmake -DUSE_CUDA ..
    make -j$(nproc)

With the `USE_CUDA` compile option activated, the available GPU architectures will be detected automatically and the code will be compiled accordingly. 

## Compiler compatibility

The standard compilation *without* CUDA should be unproblematic, but generating he CUDA-accelerated version can be more complicated. A frequent reason for difficulties is an **incompatibility of host and device compiler versions**. The compiler requirements for different CUDA distributions are summarized in a table [in this gist](https://gist.github.com/ax3l/9489132#nvcc). 


It can occur that the version of the host compiler, e.g. `g++`, is too recent for the installed CUDA version. If, for example, the output of `g++ --version` (entered on the command line prompt) yields `9.3.0` and `nvcc --version` gives `V10.2.89`, then the default `g++` compiler cannot be used. As indicated in the table referenced above, that CUDA version needs a `g++` version 8 or lower.
In such a configuration, `tetmag`'s attempt to use the CUDA compiler `nvcc` will fail and produce an error message like this:

    #error -- unsupported GNU version! gcc versions later than 8 are not supported!
To solve this problem, an older version of the host compiler must be installed. This can be done without necessarily downgrading the standard compiler, e.g., by installing `g++` version 7 in addition to the default `g++` version 9.3.0. 

Once a compatible compiler is available, `tetmag` needs to know where to find it. This information can be passed with the flag[^1] `TETMAG_HOST_COMPILER`. Assuming that the `g++-7` compiler is located in `/usr/bin/g++`, the compilation would be done with

    cd tetmag/build
    cmake -DUSE_CUDA -DTETMAG_HOST_COMPILER="/usr/bin/g++-7" ..
    make -j$(nproc)

[^1]:
    A convenient way to set parameters and options for the subsequent compilation with `cmake` and  `make` is to use [`ccmake`](https://cmake.org/cmake/help/latest/manual/ccmake.1.html). Most of the options displayed by the `ccmake` GUI refer to external libraries. The settings of the `tetmag` compilation are stored in variables named `TETMAG_*`.

## Additional software
Although `tetmag` is a standalone software that -unlike plugin-type applications- does not require a specific framework to operate, its workflow involves a series of processes that must be handled in parts by other applications. This concerns two categories of operations:

- **Preprocessing:** Definition of the sample geometry and finite-element mesh generation
- **Postprocessing:** Visualization and analysis of the results 

The following software can accomplish these tasks:

- [ParaView](https://www.paraview.org) - Data visualization 
- [gmsh](https://gmsh.info) or [netgen](https://ngsolve.org) - Finite-element mesh generation
- [FreeCAD](https://www.freecadweb.org) - Design of three-dimensional objects

For an efficient use of `tetmag`, it is recommended to install the software listed above on the machines where the simulations will be prepared and where the results will be analyzed. The required pre- and postprocessing steps and the interaction of `tetmag` with this external software will be discussed in the [tutorial](userguide.md) and in the [examples](examples.md) section. 
