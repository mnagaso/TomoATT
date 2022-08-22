# Inversion test 

This is a test setup for inversion calculation.

1. Run all cells of `make_test_model.ipynb` for creating
    - source, receiver file
    - true model
    - initial model

2. then run TOMOATT forward with `input_params_pre.yml` for calculating the true arrival times at the stations
-> this will output src_rec_result.dat file which includes the arrival time at each station

3. run TOMOATT in inversion mode with `input_params.yml`.
 