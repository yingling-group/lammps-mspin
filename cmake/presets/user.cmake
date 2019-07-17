set(PKG_ASPHERE ON CACHE BOOL "" FORCE)
set(PKG_BODY ON CACHE BOOL "" FORCE)
set(PKG_CLASS2 ON CACHE BOOL "" FORCE)
set(PKG_COLLOID ON CACHE BOOL "" FORCE)
set(PKG_DIPOLE ON CACHE BOOL "" FORCE)
set(PKG_GRANULAR ON CACHE BOOL "" FORCE)
set(PKG_KSPACE ON CACHE BOOL "" FORCE)
set(PKG_MANYBODY ON CACHE BOOL "" FORCE)
set(PKG_MC ON CACHE BOOL "" FORCE)
set(PKG_MISC ON CACHE BOOL "" FORCE)
set(PKG_MOLECULE ON CACHE BOOL "" FORCE)
set(PKG_PERI ON CACHE BOOL "" FORCE)
set(PKG_REPLICA ON CACHE BOOL "" FORCE)
set(PKG_RIGID ON CACHE BOOL "" FORCE)
set(PKG_SHOCK ON CACHE BOOL "" FORCE)
set(PKG_SNAP ON CACHE BOOL "" FORCE)
set(PKG_SRD ON CACHE BOOL "" FORCE)
set(PKG_OPT ON CACHE BOOL "" FORCE)
set(PKG_CORESHELL ON CACHE BOOL "" FORCE)
set(PKG_QEQ ON CACHE BOOL "" FORCE)

set(PKG_COMPRESS ON CACHE BOOL "" FORCE)

# set(DOWNLOAD_KIM ON CACHE BOOL "" FORCE)
# set(PKG_KIM ON CACHE BOOL "" FORCE)

# set(DOWNLOAD_LATTE ON CACHE BOOL "" FORCE)
# set(PKG_LATTE ON CACHE BOOL "" FORCE)

set(PKG_LIB ON CACHE BOOL "" FORCE)
# set(PKG_MEAM ON CACHE BOOL "" FORCE)

set(PKG_MPI ON CACHE BOOL "" FORCE)

# for spin lattice
# set(PKG_SPIN ON CACHE BOOL "" FORCE)

# set(PKG_POEMS ON CACHE BOOL "" FORCE)
set(PKG_PYTHON ON CACHE BOOL "" FORCE)
set(PKG_VOROFFOI ON CACHE BOOL "" FORCE)


set(PKG_USER OFF CACHE BOOL "" FORCE)
set(PKG_USER-BOCS OFF CACHE BOOL "" FORCE)
set(PKG_USER-CGDNA OFF CACHE BOOL "" FORCE)
set(PKG_USER-CGSDK OFF CACHE BOOL "" FORCE)
set(PKG_USER-DIFFRACTIOFF OFF CACHE BOOL "" FORCE)
set(PKG_USER-DPD OFF CACHE BOOL "" FORCE)
set(PKG_USER-DRUDE OFF CACHE BOOL "" FORCE)
set(PKG_USER-EFF OFF CACHE BOOL "" FORCE)
set(PKG_USER-FEP OFF CACHE BOOL "" FORCE)

# install libnetcdf first
set(PKG_USER-H5MD OFF CACHE BOOL "" FORCE)

set(PKG_USER-MANIFOLD OFF CACHE BOOL "" FORCE)
set(PKG_USER-MEAMC OFF CACHE BOOL "" FORCE)
set(PKG_USER-MESO OFF CACHE BOOL "" FORCE)
set(PKG_USER-MGPT OFF CACHE BOOL "" FORCE)
set(PKG_USER-MISC OFF CACHE BOOL "" FORCE)
set(PKG_USER-MOFFF OFF CACHE BOOL "" FORCE)
set(PKG_USER-MOLFILE OFF CACHE BOOL "" FORCE)

set(PKG_USER-PHONON OFF CACHE BOOL "" FORCE)
set(PKG_USER-QTB OFF CACHE BOOL "" FORCE)
set(PKG_USER-REAXC OFF CACHE BOOL "" FORCE)
set(PKG_USER-SDPD OFF CACHE BOOL "" FORCE)
set(PKG_USER-SMTBQ OFF CACHE BOOL "" FORCE)
set(PKG_USER-SPH OFF CACHE BOOL "" FORCE)
set(PKG_USER-TALLY OFF CACHE BOOL "" FORCE)
set(PKG_USER-UEF OFF CACHE BOOL "" FORCE)

# needs libraries to install
set(PKG_USER-ATC OFF CACHE BOOL "" FORCE)
set(PKG_USER-AWPMD OFF CACHE BOOL "" FORCE)
set(PKG_USER-COLVARS OFF CACHE BOOL "" FORCE)
set(PKG_USER-LB OFF CACHE BOOL "" FORCE)
set(PKG_USER-QMMM OFF CACHE BOOL "" FORCE)
set(PKG_USER-QUIP OFF CACHE BOOL "" FORCE)
set(PKG_USER-VTK OFF CACHE BOOL "" FORCE)

# needs GSL library
set(PKG_MSCG OFF CACHE BOOL "" FORCE)
set(PKG_USER-PLUMED OFF CACHE BOOL "" FORCE)

# needs netcdf library
set(NETCDF_INCLUDE_DIRS /home/akhlak/programs/miniconda3/pkgs/libnetcdf-4.6.1-h13459d8_0/include CACHE STRING "" FORCE)
set(NETCDF_LIBRARY ../../miniconda3/pkgs/libnetcdf-4.6.1-h13459d8_0/lib/libnetcdf.so.13 CACHE STRING "" FORCE)
set(PKG_USER-NETCDF OFF CACHE BOOL "" FORCE)

# needs Eigen3
set(PKG_USER-SMD OFF CACHE BOOL "" FORCE)

# needs TBB MLK
set(PKG_USER-INTEL OFF CACHE BOOL "" FORCE)

# build failed
set(PKG_MPIIO OFF CACHE BOOL "" FORCE)
set(PKG_KOKKOS OFF CACHE BOOL "" FORCE)


# turn on the MACROSPIN package we developed
set(PKG_MSPIN ON CACHE BOOL "" FORCE)

# turn on OMP version
set(BUILD_OMP ON CACHE BOOL "" FORCE)
# takes a long time to compile
set(PKG_USER-OMP ON CACHE BOOL "" FORCE)

# turn on MPI version
# make sure you installed openmpi or mpich 
# sudo apt-get install libopenmpi-dev should work
set(BUILD_MPI ON CACHE BOOL "" FORCE)

# path to installation
set(CMAKE_INSTALL_PREFIX ./ CACHE STRING "" FORCE)

# set(LAMMPS_MACHINE omp CACHE STRING "" FORCE)
set(LAMMPS_MACHINE mpi CACHE STRING "" FORCE)

# turn on GPU support
set(PKG_GPU OFF CACHE BOOL "" FORCE)
set(GPU_API cuda CACHE STRING "" FORCE)
set(GPU_ARCH sm_75 CACHE STRING "" FORCE)
# set(GPU_PREC=mixed CACHE STRING "" FORCE)

