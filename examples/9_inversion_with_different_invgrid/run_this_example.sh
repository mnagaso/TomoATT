#!/bin/bash

mkdir OUTPUT_FILES

# compute synthetic traveltime in the ckb model
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yml

# do ckb inversion using abs data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_good_two_invgrid.yml
