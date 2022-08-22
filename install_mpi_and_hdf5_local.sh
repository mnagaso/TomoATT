#!/bin/bash

# compile and install hdf5 parallel on linux (tested only on ubuntu 22.04)

# every file will be placed in external_libs
cd ./external_libs

# make a local install pass
mkdir local_mpi_hdf5

# install openmpi

# download source 
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.3.tar.gz

# Extract it
tar -xvf openmpi-4.1.3.tar.gz && cd openmpi-4.1.3

# configure

# cuda aware mpi
./configure --prefix=$(pwd)/../local_mpi_hdf5  --with-cuda
# without cuda extension
#./configure --prefix=$(pwd)/../local_mpi_hdf5



# make
make -j16 && make install

#
cd ../

# download hdf5 source
wget https://gamma.hdfgroup.org/ftp/pub/outgoing/hdf5/snapshots/v112/hdf5-1.12.2-1.tar.gz

#Extract the downloaded directory
tar -xvf hdf5-1.12.2-1.tar.gz && cd hdf5-1.12.2-1

# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=$(pwd)/../local_mpi_hdf5/bin/mpicc CXX=$(pwd)/../local_mpi_hdf5/bin/mpicxx ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

# make
make -j16 && make install

# now openmpi and hdf5 executables are in external_libs/local_mpi_hdf5/bin
