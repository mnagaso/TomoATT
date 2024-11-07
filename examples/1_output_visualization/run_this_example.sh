#!/bin/bash

mkdir 1_src_rec_files
mkdir 2_models
mkdir 3_input_params
mkdir OUTPUT_FILES


NPROC=6

# mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i 3_input_params/input_params_signal.yaml

mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i 3_input_params/input_params_inv.yaml
