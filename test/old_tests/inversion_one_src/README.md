# Inversion test 

This is a test setup for inversion calculation.

1. Run all cells of `make_test_model.ipynb` or `make_test_model.py` for creating
    - source, receiver file (src_rec_test.dat)
    - true model (test_model_true.h5)
    - initial model (test_model_init.h5)

2. then run TOMOATT forward with `input_params_pre.yml` for calculating the true arrival times at the stations
-> this will output src_rec_test_out.dat file which includes the true arrival times
```bash
mpirun --oversubscribe -np 8 ../../build/TOMOATT ./input_params_pre.yml 
```

3. run TOMOATT in inversion mode with `input_params.yml`.
-> this will output src_rec_test_out.dat file which includes the true arrival times
```bash
mpirun --oversubscribe -np 8 ../../build/TOMOATT ./input_params_pre.yml 
```

4. for visualizing the result files
```bash
paraview out_data_sim_0.xmf
```

0 is the source id. The kernel and model fields are output in the file of 0th source. 


 