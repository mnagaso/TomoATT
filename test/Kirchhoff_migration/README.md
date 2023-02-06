# Example for Kirchhoff migration mode

This example shows how to use the Kirchhoff migration mode in TOMOATT.

## compile TOMOATT with Kirchhoff migration mode

Uncomment "src/TOMOATT_km.cxx" and "src/TOMOATT_2d_precalc.cxx" in CMakeLists.txt to compile the solver_only executable as below.

```cmake
# add one by one
set(APP_SOURCES
  src/TOMOATT.cxx
  src/TOMOATT_km.cxx
  src/TOMOATT_2d_precalc.cxx
  #src/SrcRecWeight.cxx
  )
```

Then recompile the code in the build directory.

``

## run the Kirchhoff migration executable

Before running the km executable, you need to prepare the input files by running:
    
```python
python make_test_model.py
```

This creates a src_rec_file.dat which includes 8 sources without any receiver.

Then run the solver_only executable as below.

```bash
mpirun -n 8 ../../build/bin/TOMOATT_km -i input_params_pre.yml
```

Below is the output files from this example:
- OUTPUT_FILES : directory containing the output files from forward simulation
  - src_rec_file_forward.dat : new src rec file which includes calculated traveltimes at the receivers
  - out_data_sim.xmf : index file for paraview visualization.
  - out_data_sim.h5 : this includes traveltime fields from forward simulation. 

Contents of out_data_sim.h5:
```
/                        Group
/src_0                   Group 
/src_0/T_res_inv_0000    Dataset {28050} (1d travel time array for paraview visualization)
/src_0/T_res_merged_inv_0000 Dataset {10, 50, 50} (3d travel time array (already merged) for users)
```

- OUTPUT_FILES_back : directory containing the output files from backward simulation
  - out_data_sim.h5 : this includes traveltime fields from backward simulation.
  - out_data_sim.xmf : index file for paraview visualization.