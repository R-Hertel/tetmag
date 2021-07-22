# Getting started
This page describes the basic usage of `tetmag` and the content of the input files required to define the setup to be simulated.


Once the code has been successfully compiled, the executable file `tetmag` should be copied or moved to a directory that is included in `$PATH`, such as, e.g., `$HOME/bin/` or `/usr/local/bin/`, so that the programm can be run by typing `tetmag` in the command line interface.

Each simulation should have its dedicated directory. That directory must contain three specific files: `material001.dat`, `simulation.cfg`, and `<name>.msh`. The following sections explain the meaning and the content of these files.

##   **Material properties**

As the name suggests, the `material001.dat` file contains information on the material-specific properties. It is an ASCII file that can be edited according the the problem defintion. Its content is mostly self-explanatory for users who are familiar with micromagnetic theory:

````
A = 1.3e-11 # exchange anisotropy constant in J/m
Js = 1.0 # saturation polarization in T (Js = mu0 * Ms)
Ku =  500 # uniaxial anisotropy constant in J/m^3
Ks =  0  # surface anisotropy constant in J/m^2
theta_u = 90 # polar angle of uniaxial anisotropy axis in degs
phi_u = 0 # azimuthal angle of uniaxial anisotropy axis in degs
Kc1 = 0 # first cubic anisotropy constant in J/m^3
Kc2 = 0 # second cubic anisotropy constant in J/m^3
phi_Euler = 0  # cubic anisotropy axes: first Euler angle in degs
theta_Euler = 0 # .. second Euler angle in degs
psi_Euler = 0 # .. third Euler angle in degs
D = 0 # Bulk-type DMI constant in J/m^2
````

 The parameters $A$, $J_s=\mu_0M_s$, $D$, $K_u$, $K_{c1}$, and $K_{c2}$ describe the intrinsic properties of the ferromagnetic material. This input data is supplied to `tetmag` by the `material001.dat` file. Numerical values of these parameters for various magnetic materials are reported in the scientific literature. The example above shows a typical set of parameters for Permalloy (Ni$_{81}$Fe$_{19}$).

 The direction of the easy axis is defined by [two angles](https://en.wikipedia.org/wiki/Spherical_coordinate_system#/media/File:3D_Spherical.svg): $\theta_u$ is the polar angle, i.e., the angle that the easy axis encloses with the $z$ axis, and $\phi_u$ is the angle that the easy axes, when projected onto the $xy$-plane, encloses with the $x$-axis.

 Magnetic materials with cubic anisotropy have three easy axes that define an orthogonal tripod. The orientation of this orthogonal set of axes with respect to the coordinate axis system $xyz$ is described by the Euler angles $\phi_E$, $\theta_E$ and $\psi_E$. The rotation of the frame is performed according to the [$zx'z'$ convention](https://en.wikipedia.org/wiki/Euler_angles#Conventions_by_intrinsic_rotations): first a rotation by $\phi_E$ around the $x$ axis, followed by a rotation by $\theta_E$ around the $z'$ axis, and finally a rotation by $\psi_E$ around $x'$.

The strength of the Dzyaloshinskii-Moriya interaction (DMI) is described by the material parameter $D$. In helicoidal magnetic materials with intrinsic DMI, this bulk-type interaction gives rise to an energy density of the form $e_\text{DMI}=D\cdot\vec{M}\left(\nabla\times\vec{M}\right)$. Interfacial DMI, which can paly a role in ultrathin magnetic films, is currently not implemented.
 

 The current version of `tetmag` only allows for a single material type.  An extension enabling simulations of heterogeneous magnetic systems might be added in a future release.

!!! note 
    The number "001" in the file name `material001.dat` has no meaning in the current version of `tetmag`. Neverthless, the file containing the material data must have *exactly this file name*.  


##  **Simulation parameters**

 The content of the `simulation.cfg` file defines the setup of the simulation. This file  is read by `tetmag` immediately at the start. It contains information on various simulation settings: specifications about external magnetic fields, spin-polarized electrical currents, initial conditions, and technical parameters that govern the flow of the simulation. 

The short introduction on this page discusses only a small set of basic configuration settings. An exhaustive description of the options that can be specified in the `simulation.cfg` file is given in the [user guide](equations.md).

This is an example of a simple `simulation.cfg` file:

````
name = cube 
mesh type = msh

alpha = 0.1
scale = 1.e-8

initial state = random
duration = 1000 # simulation time in ps
torque limit = 1.e-4 # termination criterion: stop if maximum torque is smaller
solver type = CPU 

time step = 0.1  #  stray field freeze time in ps
    console stride = 10 # interval size in ps between outputs in console 
    log stride = 1
    config stride = 50

external field = 00.0 #  Hext in mT
	theta_H = 90 # polar angle of the field direction in degree
	phi_H = 0 # azimuthal angle
```` 
