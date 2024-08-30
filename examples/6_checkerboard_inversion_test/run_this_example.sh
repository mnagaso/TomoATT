#!/bin/bash

# --------- STEP 1, generate input parameter file--------- 
# run all cells of "1_generate_input_params.ipynb" to obtain the input_params files:
# a) input_params/input_params_signal.yaml
# b) input_params/input_params_inv_abs.yaml

mkdir -p input_params
mkdir -p OUTPUT_FILES

# --------- STEP 1.1 generate input_params_signal.yaml ---------
# generate input_params_signal.yaml
input_fname=../0_generate_files_for_TomoATT/3_input_params/0_input_params_forward_simulation.yaml
output_fname=./input_params/input_params_signal.yaml
cp $input_fname $output_fname

# update src_rec_file
pta setpar $output_fname source src_rec_file ../0_generate_files_for_TomoATT/1_src_rec_files/src_rec_config.dat

# update init_model_path
pta setpar $output_fname model init_model_path ../0_generate_files_for_TomoATT/2_models/model_ckb_N61_61_61.h5

# update_output_path
pta setpar $output_fname output_setting output_dir ./OUTPUT_FILES/OUTPUT_FILES_signal

# --------- STEP 1.2 generate input_params_inv_abs.yaml ---------

# generate input_params_inv_abs.yaml
input_fname=../0_generate_files_for_TomoATT/3_input_params/1_input_params_inversion.yaml
output_fname=./input_params/input_params_inv_abs.yaml
cp $input_fname $output_fname

# update src_rec_file
pta setpar $output_fname source src_rec_file OUTPUT_FILES/OUTPUT_FILES_signal/src_rec_file_forward_noisy.dat

# update init_model_path
pta setpar $output_fname model init_model_path ../0_generate_files_for_TomoATT/2_models/model_init_N61_61_61.h5

# update_output_path
pta setpar $output_fname output_setting output_dir OUTPUT_FILES/OUTPUT_FILES_inv_abs

# --------- STEP 2, ckb inversion --------- 
# compute synthetic traveltime in the ckb model
NPROC=8

# For Linux and Mac
# do forward simulation
# mpirun -np $NPROC ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

# do ckb inversion using abs data
# mpirun -np $NPROC ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yaml

# For WSL
# do forward simulation
mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_signal.yaml

python assign_gaussian_noise.py

# do ckb inversion using abs data
mpirun -n $NPROC --allow-run-as-root --oversubscribe ../../build/bin/TOMOATT -i input_params/input_params_inv_abs.yaml

# --------- STEP 3, plot result --------- 
# run all cells of '3_plot_ckb_model.ipynb' and '4_plot_inversion_result.ipynb' to show the results