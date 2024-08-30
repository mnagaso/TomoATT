#!/bin/bash

mkdir -p OUTPUT_FILES

# --------- STEP 2, ckb inversion --------- 
# compute synthetic traveltime in the ckb model
NPROC=8

# For Linux and Mac
# do forward simulation
# mpirun -np $NPROC ../../build/bin/TOMOATT -i 3_input_params/input_params_signal.yaml

# do ckb inversion using abs data
# mpirun -np $NPROC ../../build/bin/TOMOATT -i 3_input_params/input_params_inv_abs.yaml

# For WSL
# do forward simulation
mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i 3_input_params/input_params_signal.yaml

# do ckb inversion using abs data
# mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i 3_input_params/input_params_inv_abs.yaml

# --------- STEP 3, plot result --------- 
# run all cells of '3_plot_ckb_model.ipynb' and '4_plot_inversion_result.ipynb' to show the results