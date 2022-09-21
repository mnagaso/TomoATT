# commands for compiling and running on HPCs


## Fugaku

### 0. start interactive job on Fugaku (for accessing arch64 environment)

```

```bash 
pjsub --interact -g hp220155 -L "node=1" -L "rscgrp=int" -L "elapse=1:00:00" --sparam "wait-time=600" -x PJM_LLIO_GFSCACHE=/vol0004│17640
 --no-check-directory     
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
pjsub --interact -g hp(groupid) -L "node=1" -L "rscgrp=int" -L "elapse=1:00:00" --sparam "wait-time=600" -x PJM_LLIO_GFSCACHE=/vol0004
```

### load compiler
```bash
spack load /nphnrhl

. /vol0004/apps/oss/spack/share/spack/setup-env.sh

spack load ...
```

### 1. Download openmpi and hdf5 source code and compile it

Run the commands below (same contents in the script file)


```bash
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
#./configure --prefix=$(pwd)/../local_mpi_hdf5  --with-cuda
# without cuda extension
./configure --prefix=$(pwd)/../local_mpi_hdf5

# make
make -j16 && make install

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