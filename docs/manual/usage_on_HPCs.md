# commands for compiling and running on HPCs


## Fugaku @ RIKEN

### 0. start interactive job on Fugaku (for accessing arch64 environment)

```bash 
pjsub --interact -g hp220155 -L "node=1" -L "rscgrp=int" -L "elapse=1:00:00" --sparam "wait-time=600" -x PJM_LLIO_GFSCACHE=/vol0004 --no-check-directory     
```


### 1. Load necessary modules 
```bash
# prepare spack env
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
# load gcc 11.2.0 and Fujitsu mpi
spack load /nphnrhl /cvur4ou
```


### 2. Download hdf5 source code and compile it
```bash

# every file will be placed in external_libs
cd ./external_libs

# make a local install pass
mkdir local_mpi_hdf5

# download hdf5 source
wget https://gamma.hdfgroup.org/ftp/pub/outgoing/hdf5/snapshots/v112/hdf5-1.12.2-1.tar.gz

#Extract the downloaded directory
tar -xvf hdf5-1.12.2-1.tar.gz && cd hdf5-1.12.2-1

# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=mpicc CXX=mpic++ ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

# make
make -j16 && make install

# now hdf5 executables are in external_libs/local_mpi_hdf5/bin
```

### 3. Compile TomoATT
```bash
# cd to TomoATT directory
cd ../..

# make a build directory
mkdir build

# compile TomoATT
cd build
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5

make -j16
```

### 4. terminalte interactive job
`Ctrl + D`




## ASPIRE 1 @ NSCC

### 0. load necessary modules
```bash
module load intel/19.0.0.117 intelmpi autoconf/2.69
```

### 1. Download hdf5 source code and compile it
```bash

# on your LOCAL MACHINE (intel compilor requires hdf5 version 1.13.0 or higher)
#wget https://gamma.hdfgroup.org/ftp/pub/outgoing/hdf5/snapshots/v113/hdf5-1.13.0-rc6.tar.gz
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz


# then upload it to NSCC (for example)
scp hdf5-1.13.2.tar.gz aspire1:~/(where TomoATT placed)/external_libs/

# on ASPIURE 1
cd external_libs

mkdir local_mpi_hdf5

# extract tar file and cd to the directory
tar -xvf hdf5-1.13.2.tar.gz && cd hdf5-1.13.2

./autogen.sh

# configure the code
CC=mpiicc  CFLAGS="-I/app/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include -L/app/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib -fPIC -O3 -xHost -ip -fno-alias -align"  CXX=mpiicpc  CXXFLAGS="-I/app/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/include -L/app/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib -fPIC -O3 -xHost -ip -fno-alias -align" ./configure --with-pic --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

# make and install to the prefix
make -j16 && make install

```

### 2. Compile TomoATT
```bash
# cd to TomoATT directory
cd ../..

# make a build directory
mkdir build

# compile TomoATT
cd build
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5

make -j16
```
