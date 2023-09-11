#!/bin/bash

# python make_test_model.py

mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yml

# invert for ckb model and relocation using absolute traveltime
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs.yml

mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_reloc_abs_based_on_location.yml

# invert for ckb model and relocation using absolute traveltime
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yml

# # invert for ckb model use mode 3
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_mode3.yml
