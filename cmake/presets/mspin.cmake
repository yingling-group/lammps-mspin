# Preset that turns on just a few, frequently used packages along with the MSPIN package.
# Please note that, the RIGID package must be enabled for MSPIN to work.
# USAGE: cmake ../cmake/ -C ../cmake/presets/mspin.cmake

set(ALL_PACKAGES KSPACE MOLECULE RIGID MSPIN EXTRA-DUMP)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

# Update as necessary
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/mspin" CACHE PATH "Default install path" FORCE)
set(LAMMPS_MACHINE serial CACHE STRING "" FORCE)

# Turn on MPI support
# Make sure you installed openmpi or mpich
# sudo apt-get install libopenmpi-dev should work
set(MPI_CXX "icpx" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
set(BUILD_MPI ON CACHE BOOL "" FORCE)
set(LAMMPS_MACHINE mpi CACHE STRING "" FORCE)

# Turn on GPU support
# set(PKG_GPU ON CACHE BOOL "" FORCE)
# set(GPU_API cuda CACHE STRING "" FORCE)
# set(LAMMPS_MACHINE cuda CACHE STRING "" FORCE)

# Uncomment the gencode line corresponding to your GPU card
# Turing: RTX 2070
# set(GPU_ARCH sm_75 CACHE STRING "" FORCE)

# Maxwell: Quadro M6000 , GeForce 900, GTX-970, GTX-980, GTX Titan X
# set(GPU_ARCH sm_52 CACHE STRING "" FORCE)

# Pascal: GTX 1080, GTX 1070, GTX 1060, GTX 1050, GTX 1030 (GP108),
# GT 1010 (GP108) Titan Xp, Tesla P40, Tesla P4
# set(GPU_ARCH sm_61 CACHE STRING "" FORCE)

# Kepler: Tesla K40
# set(GPU_ARCH sm_35 CACHE STRING "" FORCE)

# (Module) Load the latest GCC, MPICH, CUDA to your terminal session.
# Set cuda lib path, if it is not in the default location.
# Or use -DCMAKE_LIBRARY_PATH=$CUDA_HOME/lib64/stubs as cmake argument
# set(CUDA_CUDA_LIBRARY /usr/local/apps/cuda/9.1/lib64/stubs/libcuda.so CACHE STRING "" FORCE)
