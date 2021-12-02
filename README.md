# LAMMPS-MSPIN
A complete LAMMPS source branch containing `MSPIN` plugin for atomistic
simulations of magnetic nanoparticles.

**Note**: If you are interested only in the source code relevant to the `MSPIN`
plugin, please see the [plugin branch](/yingling-group/lammps-mspin/tree/plugin "plugin branch").

Developed for simulation method published in:
>A.U. Mahmood and Y.G. Yingling. *All-Atom Simulation Method for Zeeman Alignment
and Dipolar Assembly of Magnetic Nanoparticles.* (2021)

Please cite paper above if used.

Author:   
Akhlak Mahmood   
Yingling Group  
North Carolina State University

# Installation Instructions

1. Clone the main branch to your local machine.
    ```
    git clone --branch main https://github.com/yingling-group/lammps-mspin.git
    ```

2. CD to the cloned directory and create a build directory.
    ```
    cd lammps-mspin
    mkdir build
    ```
3. Run CMake with the `user.cmake` preset from the build directory.
    ```
    cd build
    cmake3 -C ../cmake/presets/user.cmake ../cmake
    ```

4. Build.
    ```
    make -j4
    ```
