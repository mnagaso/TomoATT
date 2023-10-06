#!/bin/bash

# compute synthetic traveltime in the ckb model
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_signal.yml

mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params_inv_abs_cr_reloc_abs_cr.yml

