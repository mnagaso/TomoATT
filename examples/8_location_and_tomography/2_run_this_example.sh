#!/bin/bash

# --------- STEP 1, generate input parameter file--------- 
# run all cells of "1_generate_input_params.ipynb" to obtain the input_params files:
# a) input_params/input_params_signal.yaml
# b) input_params/input_params_inv_abs.yaml

# --------- STEP 2, compute synthetic traveltime in the ckb model --------- 
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

# --------- STEP 3, --------- 
# run all cells of "3_generate_obs_src_rec_data.ipynb" to obtain observation traveltime with deviated ortime and location:
# a) src_rec_obs.dat

# --------- STEP 4, update model parameters and relocate earthquake simultaneously --------- 
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_reloc_abs.yaml

# --------- STEP 5, ---------
# run all cells of '4_compare_location_result.ipynb' to compare the relocated result with the true result
# run all cells of '5_plot_location_result.ipynb' to show the relocated result

# --------- STEP 6, ---------
# run all cells of '6_plot_ckb_model.ipynb' and '7_plot_inversion_result.ipynb' to show the results

