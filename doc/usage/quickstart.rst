Getting started with tetmag
===========================

Simulations with tetmag will generate a number of output files. To avoid confusion, each simulation project should be run in a dedicated directory.
To run a simulation with ``tetmag``, three specific files are required, which must be stored in the working directory:

* ``simulation.cfg`` 

*  ``material001.dat``

* ``<name>.msh`` or ``<name>.vtk`` or ``<name>.vtu``

Here ``<name>`` is the name of your simulation project. The first two files are in human-readable ASCII format. They provide input data and information about the simulation parameters and the material properties, respectively. The third file in the list contains data on the finite-element mesh of the simulation geometry. The mesh data can be supplied in GMSH, VTK or VTU format.

The following sections will provide a few examples with a step-by-step description on how to conduct various simulations with ``tetmag``, thereby demonstrating the use and importance of these input files. The examples will also discuss how the computed output files can be analyzed.

The source code is available under AGPL-3.0 license at https://github.com/R-Hertel/tetmag .
