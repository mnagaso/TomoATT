# Eikonal equation solver based Tomography library for Adjoint-State Travel-Time tomography : TomoATT


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
