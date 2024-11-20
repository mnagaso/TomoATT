# Inversion test 

This is an example for inversion calculation with small grid.

![](img/result_fun.png)

1. Run all cells of `make_test_model.ipynb` or python script `make_test_mode.py` for creating necessary input files: 
    - source, receiver file (src_rec_test.dat)
    - true model (test_model_true.h5)
    - initial model (test_model_init.h5)

2. then run TOMOATT forward with `input_params_pre.yml` for calculating the true arrival times at the stations.
The calculated travel time at the stations is saved in the file `src_rec_test_out.dat`  
Following command will run the forward simulation with the true model
``` bash
mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params_pre.yml
```
Volumetric output data is saved in the file `OUTPUT_FILES/out_data_sim_0.h5`.
This file may be visualized by paraview with opening the index file `OUTPUT_FILES/out_data_sim_0.xmf`.
  
3. run TOMOATT in inversion mode with `input_params.yml`, by the command
``` bash
mpirun --oversubscribe -n 8 ../../build/bin/TOMOATT -i input_params.yml
```
The volumetric output data is again saved in the file `OUTPUT_FILES/out_data_sim_0.h5`.


4. The final output data is stored in OUTPUT_FILES/final_model.h5. This file may be visualized by e.g. check_3d_out.ipynb