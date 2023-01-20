# The MSPIN Package
This package contains a `fix rigid/nvt/mspin` command that updates nanoparticle
dynamics subjected to external magnetic field and mangetic dipolar interactions.

It also contains commands to compute the externel field interaction energy,
dipolar interaction energy, and interparticle distance during simulation.

See the doc page for the `fix rigid/nvt/mspin` or the `compute mspin/energy`
or `compute mspin/distance` commands for detailed usage instructions.

Use of this package requires LAMMPS to be built with the RIGID package.

There are example scripts for using commands in this package in the
examples/mspin directory.

The authors of the package is Akhlak U. Mahmood (amahmoo3 at ncsu dot edu)
and Yaroslava G. Yingling (yara_yingling at ncsu dot edu) at North Carolina
State University, USA. Contact the authors directly if you have questions.

# Installation
CMake based installation preset.

```cmake
# Enable required packages
set(ALL_PACKAGES KSPACE MOLECULE RIGID MSPIN)
foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

# Update as necessary
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/mspin" CACHE PATH "Default install path" FORCE)
set(LAMMPS_MACHINE serial CACHE STRING "" FORCE)

# Turn on MPI support
# Make sure you installed openmpi or mpich
# apt-get install libopenmpi-dev
# set(MPI_CXX "icpx" CACHE STRING "" FORCE)
# set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
# set(BUILD_MPI ON CACHE BOOL "" FORCE)
# set(LAMMPS_MACHINE mpi CACHE STRING "" FORCE)
```

# Update of the official code
List of all modifications:
```sh
$ git diff lammps/stable --name-only

cmake/CMakeLists.txt
doc/src/compute_mspin_distance.rst
doc/src/compute_mspin_energy.rst
doc/src/fix_rigid_mspin.rst
examples/mspin/README.md
examples/mspin/data.mspin
examples/mspin/in.mspin
examples/mspin/log.3Aug2022.g++.1
examples/mspin/log.3Aug2022.g++.4
src/MSPIN/README.md
src/MSPIN/compute_mspin.cpp
src/MSPIN/compute_mspin.h
src/MSPIN/compute_mspin_distance.cpp
src/MSPIN/compute_mspin_distance.h
src/MSPIN/fix_mspin_nh.cpp
src/MSPIN/fix_mspin_nh.h
src/MSPIN/fix_mspin_nvt.cpp
src/MSPIN/fix_mspin_nvt.h
src/Makefile
src/RIGID/fix_rigid.cpp
```

Other than adding the package specific files in the `src/MSPIN`, `doc/src` and `examples/mspin` directories, the following
changes are made to the official LAMMPS files.

- `src/RIGID/fix_rigid.cpp` is updated to allow additional arguments.
- `src/Makefile` is updated to add the name of the package to the *make* PACKAGE list.
- `cmake/CMakeList.txt` is updated to add the name of the pacakge to the *cmake* STANDARD_PACKAGES list.
- `.github` directory is removed.

# Recent Changes
**Jan 20, 2023**
- Updated to the latest LAMMPS stable branch.
[88c8b6ec6feac6740d140393a0d409437f637f8b](https://github.com/lammps/lammps/commit/88c8b6ec6feac6740d140393a0d409437f637f8b)

