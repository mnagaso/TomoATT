#!/bin/bash

# compile and install hdf5 parallel on linux (tested only on ubuntu 22.04)

# every file will be placed in external_libs
cd ./external_libs

# make a local install pass
mkdir local_mpi_hdf5_clang

# download hdf5 source
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz

#Extract the downloaded directory
tar -xvf hdf5-1.13.2.tar.gz
cd hdf5-1.13.2
./autogen.sh

# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
OMPI_CXX=clang++  OMPI_MPICC=clang \
CC=mpicc CXX=mpicxx \
CFLAGS="-stdlib=libc++ -rtlib=compiler-rt -unwindlib=libgcc" CXXFLAGS="-stdlib=libc++ -rtlib=compiler-rt -unwindlib=libgcc" \
./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5_clang

# make
OMPI_CXX=clang++  OMPI_MPICC=clang make -j16 && make install

# now openmpi and hdf5 executables are in external_libs/local_mpi_hdf5/bin
# cmake command would be something like this:
#OMPI_CXX=clang++ OMPI_MPICC=clang CC=clang CXX=clang++ cmake .. -DUSE_CUDA=True -DCMAKE_PREFIX_PATH=./../external_libs/local_mpi_hdf5_clang