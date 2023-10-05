#!/bin/bash

mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
# mpirun -n 3 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_signal.yml

mpirun -n 1 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_inv_abs_cs_reloc_abs_cr.yml

