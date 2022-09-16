# commands for compiling and running on HPCs


## Fugaku


### 1. Load necessary modules  
`module load mpi`  


### 2. Download hdf5 source code and compile it
```bash
cd TomoATT/external_libs

# download hdf5 source
wget https://gamma.hdfgroup.org/ftp/pub/outgoing/hdf5/snapshots/v112/hdf5-1.12.2-1.tar.gz

#Extract the downloaded directory
tar -xvf hdf5-1.12.2-1.tar.gz && cd hdf5-1.12.2-1

# configure
CC=mpicc CXX=mpicxx ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$(pwd)/../local_mpi_hdf5

# compile
make && make install
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