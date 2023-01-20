# The MSPIN Package
This package contains a `fix rigid/nvt/mspin` command that updates nanoparticle
dynamics subjected to external magnetic field and mangetic dipolar interactions.

It also contains commands to compute the externel field interaction energy,
dipolar interaction energy, and interparticle distance during simulation.

See the doc page for the `fix rigid/nvt/mspin` or the `compute mspin/energy`
or `compute mspin/distance` commands for detailed usage instructions.

Use of this package requires LAMMPS to be built with the RIGID package.

There are example scripts for using commands in this package in the
`examples/mspin` directory.

The authors of the package is Akhlak U. Mahmood (amahmoo3 at ncsu dot edu)
and Yaroslava G. Yingling (yara_yingling at ncsu dot edu) at North Carolina
State University, USA. Contact the authors directly if you have questions.

Developed for simulation method described in:
> A.U. Mahmood and Y.G. Yingling. *All-Atom Simulation Method for Zeeman Alignment
and Dipolar Assembly of Magnetic Nanoparticles.* **Journal of Chemical Theory and Computation** (2022) [doi:10.1021/acs.jctc.1c01253](https://doi.org/10.1021/acs.jctc.1c01253 "DOI").

# Installation

```sh
git clone https://github.com/yingling-group/lammps-mspin.git
cd lammps-mspin
mkdir build
cd build
cmake -C ../cmake/presets/mspin.cmake ../cmake
make -j4
```
Please update the `cmake/presets/mspin.cmake` preset file according to your machine's configuration before building.

# Usage
This package adds one `fix` and two `computes`. Please see the following doc files for usage details.
- doc/src/fix_rigid_mspin.rst
- doc/src/compute_mspin_distance.rst
- doc/src/compute_mspin_energy.rst

# Update of the official code
List of all modifications:
```sh
$ git diff lammps/stable --name-only

cmake/CMakeLists.txt
cmake/presets/mspin.cmake
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

- `src/RIGID/fix_rigid.cpp` updated to allow additional arguments.
- `src/Makefile` updated to add the name of the package to the *make* PACKAGE list.
- `cmake/CMakeList.txt` updated to add the name of the pacakge to the *cmake* STANDARD_PACKAGES list.
- `.github` directory removed.

# Recent Changes
**Jan 20, 2023**
- Updated to the latest LAMMPS stable branch.
[88c8b6ec6feac6740d140393a0d409437f637f8b](https://github.com/lammps/lammps/commit/88c8b6ec6feac6740d140393a0d409437f637f8b)
- Modified the source code to add and use a custom atom property for the `qm` magnetic charge values.
