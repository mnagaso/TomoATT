#!/bin/bash

# --------- STEP 1, generate input parameter file--------- 
# run all cells of "1_generate_input_params.ipynb" to obtain the input_params files:
# a) input_params/input_params_signal.yaml
# b) input_params/input_params_inv_abs.yaml

# --------- STEP 2, ckb inversion --------- 
# compute synthetic traveltime in the ckb model
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

# do ckb inversion using abs data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yaml

# --------- STEP 3, plot result --------- 
# run all cells of '3_plot_ckb_model.ipynb' and '4_plot_inversion_result.ipynb' to show the results