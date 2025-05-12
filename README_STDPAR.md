## About OpenFOAM STDPAR
OpenFOAM is a free, open source CFD software [released and developed by OpenCFD Ltd since 2004](http://www.openfoam.com/history/).
It has a large user base across most areas of engineering and science, from both commercial and academic organisations.
OpenFOAM has an extensive range of features to solve anything from complex fluid flows involving chemical reactions, turbulence and heat transfer, to acoustics, solid mechanics and electromagnetics.
[See documentation](http://www.openfoam.com/documentation)

OpenFOAM is professionally released every six months to include
customer sponsored developments and contributions from the community -
individual and group contributors, integrations
(eg, from FOAM-extend and OpenFOAM Foundation Ltd) as well as
[governance guided activities](https://www.openfoam.com/governance/).


## License

OpenFOAM is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.  See the file COPYING in this directory or
[http://www.gnu.org/licenses/](http://www.gnu.org/licenses), for a
description of the GNU General Public License terms under which you
may redistribute files.


## OpenFOAM Trademark

OpenCFD Ltd grants use of its OpenFOAM trademark by Third Parties on a
licence basis. ESI Group and OpenFOAM Foundation Ltd are currently
permitted to use the Name and agreed Domain Name. For information on
trademark use, please refer to the
[trademark policy guidelines][link trademark].

Please [contact OpenCFD](http://www.openfoam.com/contact) if you have
any questions about the use of the OpenFOAM trademark.

Violations of the Trademark are monitored, and will be duly prosecuted.


## Using OpenFOAM

If OpenFOAM has already been compiled on your system, simply source
the appropriate `etc/bashrc` or `etc/cshrc` file and get started.
For example, for the OpenFOAM-v2412 version:
```
source /installation/path/OpenFOAM-v2412/etc/bashrc
```

## Compiling OpenFOAM STDPAR

To compile the code ```nvc++``` compiler is needed.

The only libraries compiled with stdpar flag are:

```
src/OpenFOAM/Make/options
src/finiteVolume/Make/options
src/fvOption/Make/options
```
The flags used to compile these libraries are:

```
-DNVTX -cuda -stdpar=gpu -gpu=managed -gpu=cc90
```

## Executing laplacianFOAM with OpenFOAM and STDPAR

The folder TestLaplacianFoam contains three subfolders:
-   **heat_block_3D**: the test-case used. It consists of a cubic block made of a homogeneous material. The side walls are modeled as adiabatic. As the initial condition, a temperature of T = 40 K is applied inside a cube of size 0.2 × 0.2 × 0.2 located within the block, and T = 0 K elsewhere. To modify the initial condition, simply change the file setFieldsDict. The mesh is created using blockMesh. After generating the mesh with blockMesh, you need to run setFields to apply the initial condition;
-   **laplacianFoam**: the laplacianFOAM application used on GH200 machine, please take cares of the compilation flags on the *options* file: ```-DNVTX -stdpar=gpu -gpu=cc90 -gpu=managed ```
-   **log**: two log files obtained with the native OpenFOAM implementation log_2MPI_CPU and the STDPAR version log_2MPI_GPU;