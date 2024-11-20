# Compile HDF5 library with parallel option from the source

1. run the install script `./install_mpi_and_hdf5_local.sh`, which compiles and creates openmpi and hdf5 executables in `./external_libs/local_mpi_hdf5/bin`

2. compile TOMOATT by 
``` bash
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../external_libs/local_mpi_hdf5
make -j16
```

This creates TOMOATT executable at ./build/TOMOATT