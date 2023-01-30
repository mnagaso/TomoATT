# Example for Solver only mode

This example shows how to use the Solver only mode of the tool.

## compile TOMOATT with solver_only executable

Uncomment "src/TOMOATT_solver_only.cxx" in CMakeLists.txt to compile the solver_only executable as below.

```cmake
# add one by one
set(APP_SOURCES
  src/TOMOATT.cxx
  src/TOMOATT_solver_only.cxx
  #src/SrcRecWeight.cxx
  )
```

Then recompile the code in the build directory.

## run the solver_only executable

Before running the solver_only executable, you need to prepare the input files by running:
    
```python
python make_test_model.py
```

This creates a src_only_test.dat which includes 8 sources without any receiver.

Then run the solver_only executable as below.

```bash
mpirun -n 8 ../../build/bin/TOMOATT_solver_only -i src_only_test.dat
```

## check the output files
The result file can be visualize with paraview:
```bash
paraview OUTPUT_FILES/out_data_sim.xmf
```
or directly open with python etc. Please refer the check_3d_out.ipynb for more details.