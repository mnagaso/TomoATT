# Inversion test 

This is a test setup for inversion calculation with teleseismic events.

1. Run all cells of `make_test_model.ipynb` for creating
    - source, receiver file (src_rec_test.dat)
    - true model (test_model_true.h5)
    - initial model (test_model_init.h5)

2. then run TOMOATT forward with `input_params_pre.yml` for calculating the true arrival times at the stations
At first, TOMOATT calculates 2d traveltime field for each teleseismic event.
TOMOATT will find the event source is teleseismic if the event origin is outside of the simulation domain.  
Calculated 2d traveltime field is saved in OUTPUT_FILES directory with 2d_travel_time_field_{src_id}.h5  
2D solver run will be skipped if this file is found by TOMOATT.
Result travel time at the stations is saved in the file `src_rec_test_out.dat`  
Following command will run the forward simulation with the true model
``` bash
mpirun -n 8 ../../build/TOMOATT -i input_params_pre.yml
```
Volumetric output data is saved in the file `OUTPUT_FILES/out_data_sim_0.h5`.
This file may be visualized by paraview with opening the index file `OUTPUT_FILES/out_data_sim_0.xmf`.

3. run TOMOATT in inversion mode with `input_params.yml`, by the command
``` bash
mpirun -n 8 ../../build/TOMOATT -i input_params.yml
```
The volumetric output data is saved in the file `OUTPUT_FILES/out_data_sim_0.h5`.


4. For retrieving the output data from *.h5, please use the python script `utils/tomoatt_data_retrieval.py`