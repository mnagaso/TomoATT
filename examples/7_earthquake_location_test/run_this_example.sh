#!/bin/bash

mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yml

# run generate_syn_obs_src_rec_data.ipynb to deviate earthquake location and origin time
# python generate_syn_obs_src_rec_data.py

# run location using absolute traveltime data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs.yml

# run location using common receiver differential traveltime data
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_reloc_cr.yml

# run location using both data
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_reloc_abs_cr.yml

# python compare_location_result.py


