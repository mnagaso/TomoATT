# commands for compiling and running on HPCs

## Gekko @ NTU

### 1. Load necessary modules and select GNU compilers
```bash
module purge && module load cmake gnu/gcc-9.3.0
```

### 2. compiler openmpi and HDF5 with parallel option
`./install_mpi_and_hdf5_local.sh`  
will create openmpi and hdf5 executables in external_libs/local_mpi_hdf5/bin

### 3. Compile TomoATT
```bash
# make a build directory
mkdir build

# compile TomoATT
cd build
CC=/usr/local/gcc-9.3.0/bin/gcc CXX=/usr/local/gcc-9.3.0/bin/g++ cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5

make -j16
```

Then the executable TOMOATT is created in the build directory.


## Fugaku @ RIKEN


### 0. start interactive job on Fugaku (for accessing arch64 environment)

```bash 
pjsub --interact -g hp220155 -L "node=1" -L "rscgrp=int" -L "elapse=1:00:00" --mpi "max-proc-per-node=12" --sparam "wait-time=600" -x PJM_LLIO_GFSCACHE=/vol0004 --no-check-directory     
```


### 1. Load necessary modules 
```bash
# prepare spack env
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
# load Fujitsu mpi
spack load /jfzaut5
# or load gnu 11.2.0
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
CC=mpifcc CFLAGS="-Nclang" CXX=mpiFCC CXXFLAGS="-Nclang" ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

# make
make -j12 && make install

# now hdf5 executables are in external_libs/local_mpi_hdf5/bin
```

or with gnu 11.2.0
```bash

# download hdf5 source
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz
#Extract the downloaded directory
tar -xvf hdf5-1.13.2.tar.gz
cd hdf5-1.13.2
# Configure the code. (the pathes to mpicc, mpicxx should vary on the environment)
CC=mpicc CXX=mpic++ ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5
# make
make -j12 && make install

```


### 3. Compile TomoATT
```bash
# cd to TomoATT directory
cd ../..

# make a build directory
mkdir build

# compile TomoATT
cd build

CC=mpifcc CXX=mpiFCC cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5
# or for gnu 11.2.0
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5

make -j12
```

### 4. terminalte interactive job
`Ctrl + D`





## ASPIRE 1 @ NSCC

### 0. load necessary modules
```bash
module purge
export PATH=/app/gcc/9.5.0/bin:$PATH && export LD_LIBRARY_PATH=/app/gcc/9.5.0/lib:$LD_LIBRARY_PATH
module load intel/19.0.0.117
```

### 1. Download hdf5 source code and compile it
```bash

# download openmpi and hdf5 source code on your LOCAL MACHINE
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.3/src/hdf5-1.13.3.tar.gz

# then upload it to NSCC (for example)
scp -rC hdf5-1.13.3.tar.gz aspire1:~/(where TomoATT placed)/external_libs/

# on ASPIURE 1
cd external_libs

mkdir local_mpi_hdf5

# extract tar file and cd to the directory
tar -xvf hdf5-1.13.3.tar.gz && cd hdf5-1.13.3

# configure the code
CC=mpiicc  CXX=mpiicpc \
./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

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
CC=icc CXX=icpc cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5

make -j16
```
