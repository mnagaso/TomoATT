#!/bin/bash

# compile and install hdf5 parallel on linux (tested only on ubuntu 22.04)

# every file will be placed in external_libs
cd ./external_libs

# make a local install pass
mkdir local_mpi_hdf5_intel

# install openmpi

# download source
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.3.tar.gz

# Extract it
tar -xvf openmpi-4.1.3.tar.gz && cd openmpi-4.1.3

# configure

# cuda aware mpi
CC=icx CXX=icpx ./configure --prefix=$(pwd)/../local_mpi_hdf5_intel --with-cuda --disable-mpi-fortran
# without cuda extension
#./configure --prefix=$(pwd)/../local_mpi_hdf5_intel



# make
make -j16 && make install

#
cd ../

# download hdf5 source
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz

#Extract the downloaded directory
tar -xvf hdf5-1.13.2.tar.gz && cd hdf5-1.13.2

# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=$(pwd)/../local_mpi_hdf5_intel/bin/mpicc CXX=$(pwd)/../local_mpi_hdf5_intel/bin/mpicxx ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5_intel

# make
make -j16 && make install

# now openmpi and hdf5 executables are in external_libs/local_mpi_hdf5/bin
