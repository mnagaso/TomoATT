# Inversion test 

This is a test setup for inversion calculation.

1. Run all cells of `make_test_model.ipynb` for creating
    - source, receiver file
    - true model
    - initial model

2. then run TOMOATT forward with `input_params_pre.yml` for calculating the true arrival times at the stations
-> this will output src_rec_result.dat file which includes the arrival time at each station

3. run TOMOATT in inversion mode with `input_params.yml`.

4. at first, TOMOATT calculates 2d traveltime field for each teleseismic event.
TOMOATT will find the event source is teleseismic if the event origin is outside of the simulation domain.
Calculated 2d traveltime field is saved in OUTPUT_FILES directory with 2d_travel_time_field_{src_id}.h5
2D solver run will be skipped if this file is found by TOMOATT.
 