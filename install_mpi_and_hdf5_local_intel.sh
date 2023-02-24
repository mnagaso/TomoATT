#!/bin/bash

# compile and install hdf5 parallel on linux (tested only on ubuntu 22.04)

# every file will be placed in external_libs
cd ./external_libs

## make a local install pass
mkdir local_mpi_hdf5_intel

# download hdf5 source
#wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.3/src/hdf5-1.13.3.tar.gz
##Extract the downloaded directory
#tar -xvf hdf5-1.13.3.tar.gz
cd hdf5-1.13.3

# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=/opt/intel/oneapi/mpi/latest/bin/mpiicc CXX=/opt/intel/oneapi/mpi/latest/bin/mpiicpc \
    CFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align -I/opt/intel/oneapi/mpi/latest/include -L/opt/intel/oneapi/mpi/latest/lib:/opt/intel/oneapi/mpi/latest/lib/release" \
    CXXFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align -I/opt/intel/oneapi/mpi/latest/include -L/opt/intel/oneapi/mpi/latest/lib:/opt/intel/oneapi/mpi/latest/lib/release" \
    ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --disable-fortran --with-pic --prefix=$(pwd)/../local_mpi_hdf5_intel

# make
make -j16 && make install

# now openmpi and hdf5 executables are in external_libs/local_mpi_hdf5/bin
