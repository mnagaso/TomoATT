# Inversion test 

This is a test setup for inversion calculation.

1. `python make_test_model.py` to generate the test model
2. `mpirun --oversubscribe -n 8 ../../build/TOMOATT  -i ./input_params_pre.yml` to calculate the true traveltime at the receiver locations (output src_rec_test_out.dat)
3. `python modify_source_location.py` to modify the source location (output src_rec_test_out_modified.dat)
4. `mpirun --oversubscribe -n 8 ../../build/TOMOATT  -i ./input_params.yml` (output src_rec_test_out_modified_out.dat)
5. `python plot_check_relocated_source.py` to plot the relocated source location