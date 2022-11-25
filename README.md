# TomoATT 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)


![logo](docs/logo/TomoATT_logo_2.png)

TomoATT is a library which implements an eikonal equation solver based Adjoint-state Travel-Time Tomography for a very large-scale computation.
The details of this method is described in the published works [Ping Tong (2021)](https://doi.org/10.1029/2021JB021818) and [Jing Chen (2022)](add_here_when_published).
In this release, the slowness field is computed in spherical coordinate system. The inversion for anisotropy fields will be also available in the up comming 

This library is developped to be used for modeling a very-large domain on HPCs. For this purpose, we applied 3-layer parallelization, which are:
- layer 1: simulutaneous run parallelization (multiple seismic sources may be calculated simultaneously)
- layer 2: domain decomposition (If the number of computational nodes requires too large memory, the global domain can be separated into sub-domains and calculation fo r each sub-domain is run in a separated compute node)
- layer 3: sweeping parallelization (in each sub-domain, sweeping layers are also parallelized)

The further details of the parallelization method applied in this library are described in the paper [Miles Detrixhe and Frédéric Gibou (2016)](https://doi.org/10.1016/j.jcp.2016.06.023).

Regional events (sources within the global domain).


## dependency
- MPI v3.0 or higher  

optinal:
- HDF5 compiled with parallel option
- h5py (used in pre/post processes examples)

## to clone the repository
```
``` bash
git clone --recursive https://github.com/mnagaso/TomoATT.git
```

## to compile the code 
``` bash
mkdir build && cd build
cmake .. && make -j 8
```

or compile with cuda support
``` bash
cmake .. -DUSE_CUDA=True && make -j 8
``` 

Please check the [user manual](./docs/manual/index.md) and `examples` directory for the details.


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


### HDF5 problem
In the case you will get the error message below:
``` text
/Users/..../Codes/TomoATT/src/io.cpp:1113:5: error: use of undeclared identifier 'H5Pset_fapl_mpio'; did you mean 'H5Pset_fapl_core'?
    H5Pset_fapl_mpio(plist_id, inter_sub_comm, MPI_INFO_NULL);
    ^~~~~~~~~~~~~~~~
    H5Pset_fapl_core
```

Your HDF5 linked to TomoATT is not compiled with parallel option. Please check the HDF5 installation. If the HDF5 is properly installed, `h5pcc` command should be available in your environment.