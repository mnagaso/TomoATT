#!/bin/bash

# python make_test_model.py
# run make_test_model.ipynb


mkdir OUTPUT_FILES

# # compute synthetic traveltime in the ckb model
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_swap_signal.yml

# do ckb inversion using abs data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_swap_inv_abs.yml

# run model_visualization.ipynb to plot ckb model and inv model

# python model_visualization.py ckb
# python model_visualization.py inv