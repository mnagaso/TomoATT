#!/bin/bash

# --------- STEP 1, generate model --------- 
# run all cells of "1_make_ckb_model.ipynb" to obtain the input_params files:
# a) models/model_init_N81_81_81.h5
# b) models/model_ckb_N81_81_81.h5

# --------- STEP 2, compute synthetic traveltime in the ckb model --------- 
# mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

# link the 2D traveltime field to the OUTPUT path
# mkdir OUTPUT_FILES/OUTPUT_FILES_tele_inv
# ln -s ../OUTPUT_FILES_signal/2D_TRAVEL_TIME_FIELD OUTPUT_FILES/OUTPUT_FILES_tele_inv

# # do ckb inversion using abs data
mpirun -n 8 --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_tele.yaml

