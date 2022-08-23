# TomoATT 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)


TomoATT is a library which implements an eikonal equation solver based Adjoint-state Travel-Time Tomography for a very large-scale computation, following a published article [Ping Tong (2021)](https://doi.org/10.1029/2021JB021818) and [Jing Chen (2022)](add_here_when_published).
The slowness field and anisotropy field are computed in spherical coordinate system.

Thanks to the efficiency of an eikonal equation solver, the computation of the travel-time is very fast and requires less amount of computational resources.
As an input data for TomoATT is travel times at seismic stations, we can easily prepare a great amount of data for the computation.

This library aims to be used for modeling a very-large domain. For this purpose, 3-layer parallelization is applied, which are:
- layer 1: simulutaneous run parallelization (travel times for multiple seismic sources may be calculated simultaneously)
- layer 2: subdomain decomposition (If the number of computational nodes requires too large memory, we can separate the domain into subdomains and run each subdomain in a separate compute node)
- layer 3: sweeping parallelization (in each subdomain, sweeping layers are also parallelized)

The details of the parallelization method are described in the paper [Miles Detrixhe and Frédéric Gibou (2016)](https://doi.org/10.1016/j.jcp.2016.06.023).


## dependency
- MPI v3.0 or higher  

below is optinal:
- HDF5 (parallel IO needs to be enabled)
- h5py (used in pre/post processes examples)

## to clone
``` bash
git clone --recursive https://github.com/mnagaso/TomoATT.git
```

## to compile
``` bash
mkdir build && cd build
cmake .. && make -j 8
```

## to run an example
in build directory,
``` bash
mpirun -n 4 ./TOMOATT -i ./input_params.yml
```
Please check the `examples` directory for the details.


## FAQs.
### git submodule problem
In the case you will get the error message below:
``` text
-- /(path to your tomoatt)/TomoATT/external_libs/yaml-cpp/.git does not exist. Initializing yaml-cpp submodule ...
fatal: not a git repository (or any of the parent directories): .git
CMake Error at external_libs/CMakeLists.txt:9 (message):
  /usr/bin/git submodule update --init dependencies/yaml-cpp failed with exit
  code 128, please checkout submodules
Call Stack (most recent call first):
  external_libs/CMakeLists.txt:13 (initialize_submodule)
```

You will need to update the submodule manuary by the command:
``` bash
git submodule update --init --recursive
```

In the case that `git submodule` command doesn't work in your environment, you need to download yaml-cpp library from [its repository](https://github.com/jbeder/yaml-cpp), and place it in the `external_libs` directory,
so that this directory is placed as `external_libs/yaml-cpp`.
