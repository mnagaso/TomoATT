#!/bin/bash

# python make_test_model.py

mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yml

# run generate_syn_obs_src_rec_data.ipynb to deviate earthquake location and origin time

# invert for ckb model and relocation using absolute traveltime
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs.yml

mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs_cs_cr_reloc_abs_cr_based_on_location.yml

