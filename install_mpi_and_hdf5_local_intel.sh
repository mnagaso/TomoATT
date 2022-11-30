#!/bin/bash

# compile and install hdf5 parallel on linux (tested only on ubuntu 22.04)

# every file will be placed in external_libs
cd ./external_libs

## make a local install pass
#mkdir local_mpi_hdf5_intel
#
## install openmpi
#
## download source
#wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.3.tar.gz
#
## Extract it
#tar -xvf openmpi-4.1.3.tar.gz && cd openmpi-4.1.3
#
## configure
#
## cuda aware mpi
#CC=icx CXX=icpx ./configure --prefix=$(pwd)/../local_mpi_hdf5_intel --with-cuda --disable-mpi-fortran
## without cuda extension
##./configure --prefix=$(pwd)/../local_mpi_hdf5_intel
#
#
#
## make
#make -j16 && make install
#
##
#cd ../

## 1) szip 2.1
#export CC=icc
#export CXX=icpc
#export FC=ifort
#export CFLAGS='-O3 -xHost -ip'
#export CXXFLAGS='-O3 -xHost -ip'
#export FCFLAGS='-O3 -xHost -ip'
#wget http://sources.buildroot.net/szip/szip-2.1.tar.gz
#tar -zxvf szip-2.1.tar.gz
#cd szip-2.1
#./configure --prefix=/usr/local/szip-2.1
#make
##make check
##make install
#cd ../
#
## 2) zlib 1.2.7
#export CC=icc
#export CFLAGS='-O3 -xHost -ip'
#wget https://src.fedoraproject.org/repo/pkgs/mingw-zlib/zlib-1.2.7.tar.gz/60df6a37c56e7c1366cca812414f7b85/zlib-1.2.7.tar.gz
#tar -zxvf zlib-1.2.7.tar.gz
#cd zlib-1.2.7
#./configure --prefix=/usr/local/zlib-1.2.7
#make
##make check
##make install
#
#cd ../

#--with-szip=$(pwd)/../szip-2.1 --with-zlib=$(pwd)/../zlib-1.2.7

# download hdf5 source
#wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.3/src/hdf5-1.13.3.tar.gz
##Extract the downloaded directory
#tar -xvf hdf5-1.13.3.tar.gz
cd hdf5-1.13.3

./autogen.sh
# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=/opt/intel/oneapi/mpi/latest/bin/mpiicc CXX=/opt/intel/oneapi/mpi/latest/bin/mpiicpc \
    CFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align -I/opt/intel/oneapi/mpi/latest/include -L/opt/intel/oneapi/mpi/latest/lib:/opt/intel/oneapi/mpi/latest/lib/release" \
    CXXFLAGS="-fPIC -O3 -xHost -ip -fno-alias -align -I/opt/intel/oneapi/mpi/latest/include -L/opt/intel/oneapi/mpi/latest/lib:/opt/intel/oneapi/mpi/latest/lib/release" \
    ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --disable-fortran --with-pic --prefix=$(pwd)/../local_mpi_hdf5_intel

# make
make -j16 && make install

# now openmpi and hdf5 executables are in external_libs/local_mpi_hdf5/bin
